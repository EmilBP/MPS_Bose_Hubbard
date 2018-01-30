#include "itensor/all.h"
#include "boson.h"
#include "TimeEvolve.hpp"
#include "correlations.hpp"
#include <time.h>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>

#include "gnuplot-iostream.h"

using namespace itensor;
using std::vector;

int main() {
  int N = 15;
  int Npart = 15;
  int HilbertDim = 12;
  Real tstep = 1e-3; //time step (smaller is generally more accurate)
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


  // auto ampoSF = AutoMPO(sites);
  // for(int i = 1; i < N; ++i) {
  //   ampoSF += -1.0,"A",i,"Adag",i+1;
  //   ampoSF += -1.0,"Adag",i,"A",i+1;
  // }
  // auto H = IQMPO(ampoSF);
  //
  // auto sweeps = Sweeps(5);
  // sweeps.maxm() = 20,50,75,100,200,200;
  // sweeps.cutoff() = 1E-9;
  // sweeps.niter() = 2;
  // sweeps.noise() = 1E-7,1E-8,0;
  // println(sweeps);
  //
  // auto energy = dmrg(psi,H,sweeps,"Quiet");
  // auto psi0 = psi;


  //
  //  SETUP GATES FOR TROTTER STEP
  //
  // using Gate = BondGate<IQTensor>;
  // auto gatelist = vector<Gate>();
  //
  // for(int i = 1; i < N; ++i)
  //     {
  //     auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
  //     hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);
  //
  //     auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
  //     gatelist.push_back(g);
  //     }
  // for(int i = N-1; i >= 1; --i)
  //     {
  //     auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
  //     hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);
  //
  //     auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
  //     gatelist.push_back(g);
  //     }


  //
  //  TIME EVOLUTION USING TENSOR EXPANSION -- ORDER 10
  //
  using tensorlist = std::vector< std::vector<IQTensor> >;
  clock_t begin1  = clock();
  tensorlist expUtensor;
  for (size_t j = 0; j < ttotal/tstep; j++) {
    auto tmplist = std::vector<IQTensor>();
    double U = (1.0+4.0*tstep*j);
    // double U = 1.0;
    for (int i = 1; i <= N; ++i) {
      auto hterm = -0.5*U*tstep*Cplx_i*sites.op("N(N-1)",i);
      auto Id    = sites.op("Id",i);

      auto term = hterm;
      hterm.mapprime(1,2);
      hterm.mapprime(0,1);
      IQTensor gate_;

      for(int ord = 10; ord >= 1; --ord)
          {
          term /= ord;
          gate_ = Id + term;
          term = gate_ * hterm;
          term.mapprime(2,1);
          }

      tmplist.push_back(gate_);
    }
    expUtensor.push_back(tmplist);
  }
  clock_t end1  = clock();
  // clock_t begin2= clock();
  // auto psi1     = psi0;
  // auto psi_t1   = TimeEvolve(gatelist,expUtensor,ttotal,tstep,psi1,{"Cutoff=",cutoff,"Verbose=",true});
  // clock_t end2  = clock();
  std::cout << "Building tensor expansion = " <<  double(end1 - begin1) / CLOCKS_PER_SEC << "\n";
  // std::cout << "Time Evolution for TEBD + tensor = " <<  double(end2 - begin2) / CLOCKS_PER_SEC << "\n";


  //
  //  TIME EVOLUTION USING CUSTOM BUILT TENSORS
  //
  begin1  = clock();
  tensorlist expUtensor2;


  std::vector< std::vector<IQIndexVal> > indices;
  std::vector< std::vector<IQIndexVal> > indicesP;

  auto Tlist = std::vector<IQTensor>();
  int HD;

  for (int k = 1; k <= N; ++k) {
    auto s  = sites.si(k);
    auto sP = prime(s);

    HD = s.nblock();

    std::vector<IQIndexVal> itmp(HD);
    std::vector<IQIndexVal> itmpP(HD);

    for (size_t i = 0; i < HD; i++) {
      itmp.at(i) = s(i+1);
      itmpP.at(i) = sP(i+1);
    }

    IQTensor T(dag(s),sP);
    Tlist.push_back(T);
    indices.push_back(itmp);
    indicesP.push_back(itmpP);
  }

  for (size_t j = 0; j < ttotal/tstep; j++) {
    double U  = (1.0+4.0*tstep*j);

    int n = 0;
    for (auto& T : Tlist){

      for (size_t i = 0; i < HD; i++) {
        T.set(indices[n][i],indicesP[n][i], std::exp( -0.5*U*tstep*Cplx_i*i*(i-1) ) );
      }
      n++;
    }

    expUtensor2.push_back(Tlist);
  }


  end1          = clock();
  // begin2        = clock();
  // auto psi2     = psi0;
  // auto psi_t2   = TimeEvolve(gatelist,expUtensor2,ttotal,tstep,psi2,{"Cutoff=",cutoff,"Verbose=",true});
  // end2          = clock();
  std::cout << "Building custom tensor = " <<  double(end1 - begin1) / CLOCKS_PER_SEC << "\n";
  // std::cout << "Time Evolution for TEBD + tensor = " <<  double(end2 - begin2) / CLOCKS_PER_SEC << "\n";


  begin1 = clock();
  tensorlist expUtensor3;
  for (size_t j = 0; j < ttotal/tstep; j++) {
    double U  = (1.0+4.0*tstep*j);
    auto tmplist = std::vector<IQTensor>();


    for (int k = 1; k <= N; ++k) {
      auto s  = sites.si(k);
      auto sP = prime(s);

      std::vector<IQIndexVal> indices(HilbertDim+1);
      std::vector<IQIndexVal> indicesP(HilbertDim+1);
      for (size_t i = 0; i <= HilbertDim; i++) {
        indices.at(i) = s(i+1);
        indicesP.at(i) = sP(i+1);
      }

      IQTensor T(dag(s),sP);

      for (size_t i = 0; i <= HilbertDim; i++) {
          T.set(indices.at(i),indicesP.at(i), std::exp( -0.5*U*tstep*Cplx_i*i*(i-1) ) );
      }

      tmplist.push_back(T);
    }

    expUtensor3.push_back(tmplist);
  }
  end1          = clock();
  std::cout << "Building custom tensor = " <<  double(end1 - begin1) / CLOCKS_PER_SEC << "\n";

  //
  //  POST ANALYSIS
  //
  // std::vector<std::pair<double, double> > CF1;
  // double tcurrent = 0;
  // for (auto& psit : psi_t1){
  //   auto rho = correlationTerm(sites,psit,"Adag","A");
  //   auto fc = rho/Npart;
  //   CF1.push_back(std::make_pair(tcurrent, fc));
  //   tcurrent += tstep;
  // }
  //
  // std::vector<std::pair<double, double> > CF2;
  // tcurrent = 0;
  // for (auto& psit : psi_t2){
  //   auto rho = correlationTerm(sites,psit,"Adag","A");
  //   auto fc = rho/Npart;
  //   CF2.push_back(std::make_pair(tcurrent, fc));
  //   tcurrent += tstep;
  // }


    //
    //
    // //
    // //  SAVE RESULTS
    // //
    // std::fstream myfile;
    // std::string name = "TEalgoSFdt" + std::to_string(tstep) + "cutoff9.txt";
    // myfile.open(name,std::fstream::out);
    // for (size_t i = 0; i < CF1.size(); i++) {
    //   myfile << CF1[i].first << "\t";
    //   myfile << CF1[i].second << "\t";
    //   myfile << CF2[i].second << "\t";
    //   myfile << CF3[i].second << "\t";
    //   myfile << CF4[i].second << "\t";
    //   myfile << CF5[i].second << "\n";
    // }
    // myfile << T1 << "\t" << T2 << "\t" << T3 << "\t" << T4 << "\t" << T5 << "\n";
    // myfile.close();
    //
    //
    //  PLOT RESUTLS
    //
    // Gnuplot gp;
    // // gp << "set xrange [0:1]\n";
    // gp << "set xlabel 't'\n";
    // gp << "set ylabel 'Condensate Fraction'\n";
    //
  	// gp << "plot"
    //    << gp.file1d(CF1) << "with lines title 'Tensor 10',"
    //    << gp.file1d(CF2) << "with lines title 'Custom',"
    //    << std::endl;

  return 0;
}
