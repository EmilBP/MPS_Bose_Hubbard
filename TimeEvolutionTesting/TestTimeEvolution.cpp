#include "itensor/all.h"
#include "boson.h"
#include "TimeEvolve.hpp"
#include "correlations.hpp"
#include <time.h>

#include "gnuplot-iostream.h"

using namespace itensor;
using std::vector;

int main() {
  int N = 6;
  int Npart = 6;
  int HilbertDim = 5;
  Real tstep = 1e-2; //time step (smaller is generally more accurate)
  Real ttotal = 1.0; //total time to evolve
  Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form

  auto sites = Boson(N,HilbertDim);
  auto state = InitState(sites);
  int p = Npart;
  for(int i = N; i >= 1; --i) {
      if (p >= 1) { state.set(i,"Occ1"); p -= 1; }
      else { state.set(i,"Emp"); }
  }
  auto psi = IQMPS(state);


  auto ampoSF = AutoMPO(sites);
  for(int i = 1; i < N; ++i) {
    ampoSF += -1.0,"A",i,"Adag",i+1;
    ampoSF += -1.0,"Adag",i,"A",i+1;
  }
  auto H = IQMPO(ampoSF);

  auto sweeps = Sweeps(5);
  sweeps.maxm() = 20,50,75,100,200,200;
  sweeps.cutoff() = 1E-9;
  sweeps.niter() = 2;
  sweeps.noise() = 1E-7,1E-8,0;
  println(sweeps);

  auto energy = dmrg(psi,H,sweeps,"Quiet");
  auto psi0 = psi;

  //
  // SETUP REGULAR TIME EVOLUTION
  //
  std::vector<AutoMPO> ampolist;
  for (size_t j = 0; j < ttotal/tstep; j++) {
    auto ampo = AutoMPO(sites);
    for(int i = 1; i < N; ++i) {
      ampo += -1.0,"A",i,"Adag",i+1;
      ampo += -1.0,"Adag",i,"A",i+1;
    }
    for (int i = 0; i <= N; i++) {
      // double U = (1.0+4.0*tstep*j);
      double U = 1.0;
      ampo += U,"N(N-1)",i;
    }
    ampolist.push_back(ampo);
  }


  //
  //  SETUP GATES FOR TROTTER STEP
  //
  using Gate = BondGate<IQTensor>;
  auto gatelist = vector<Gate>();

  for(int i = 1; i < N; ++i)
      {
      auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);

      auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      gatelist.push_back(g);
      }
  for(int i = N-1; i >= 1; --i)
      {
      auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);

      auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      gatelist.push_back(g);
      }

  //
  //  REGULAR TIME EVOLUTION
  //
  auto psi1     = psi0;
  auto args     = Args("Cutoff=",1E-8,"Maxm=",300);
  clock_t begin = clock();
  auto psi_t1   = TimeEvolve(psi1,ampolist,Cplx_i*tstep,args);
  clock_t end   = clock();
  std::cout << "Runtime for regular time evolution = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";


  //
  //  TIME EVOLUTION USING TENSOR EXPANSION
  //
  auto expUtensor = vector<IQTensor>();
  for (int i = 1; i <= N; ++i) {
    auto hterm = -1.0*tstep*Cplx_i*sites.op("N(N-1)",i);
    auto Id    = sites.op("Id",i);

    auto term = hterm;
    hterm.mapprime(1,2);
    hterm.mapprime(0,1);
    IQTensor gate_;

    for(int ord = 100; ord >= 1; --ord)
        {
        term /= ord;
        gate_ = Id + term;
        term = gate_ * hterm;
        term.mapprime(2,1);
        }
    expUtensor.push_back(gate_);
  }
  auto psi2     = psi0;
  begin         = clock();
  auto psi_t2   = TimeEvolveMixed2(gatelist,expUtensor,ttotal,tstep,psi2,{"Cutoff=",cutoff,"Verbose=",true});
  end           = clock();
  std::cout << "Runtime for TEBD + tensor = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";


  //
  //  TIME EVOLUTION USING MPO
  //
  auto ampo = AutoMPO(sites);
  for(int i = 1; i <= N; ++i) {
    // double U = (1.0+4.0*tstep*j);
    double U = 1;
    ampo += U,"N(N-1)",i;
  }
  auto expUmpo = toExpH<IQTensor>(ampo,tstep*Cplx_i);

  auto psi3     = psi0;
  begin         = clock();
  auto psi_t3   = TimeEvolveMixed3(gatelist,expUmpo,ttotal,tstep,psi3,{"Cutoff=",cutoff,"Verbose=",true});
  end           = clock();
  std::cout << "Runtime for TEBD + MPO = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";


  //
  //  POST ANALYSIS
  //
  std::vector<std::pair<double, double> > CF1;
  double tcurrent = 0;
  for (auto& psit : psi_t1){
    auto rho = correlationTerm(sites,psit,"Adag","A");
    auto fc = rho/Npart;
    CF1.push_back(std::make_pair(tcurrent, fc));
    tcurrent += tstep;
  }

  std::vector<std::pair<double, double> > CF2;
  tcurrent = 0;
  for (auto& psit : psi_t2){
    auto rho = correlationTerm(sites,psit,"Adag","A");
    auto fc = rho/Npart;
    CF2.push_back(std::make_pair(tcurrent, fc));
    tcurrent += tstep;
  }

  std::vector<std::pair<double, double> > CF3;
  tcurrent = 0;
  for (auto& psit : psi_t3){
    auto rho = correlationTerm(sites,psit,"Adag","A");
    auto fc = rho/Npart;
    CF3.push_back(std::make_pair(tcurrent, fc));
    tcurrent += tstep;
  }


  //
  //  PLOT RESUTLS
  //
  Gnuplot gp;
  // gp << "set xrange [0:1]\n";
  gp << "set xlabel 't'\n";
  gp << "set ylabel 'Condensate Fraction'\n";

	gp << "plot"
     << gp.file1d(CF1) << "with lines title 'Regular',"
     << gp.file1d(CF2) << "with lines title 'Tensor',"
     << gp.file1d(CF3) << "with lines title 'MPO',"
     << std::endl;

  return 0;
}
