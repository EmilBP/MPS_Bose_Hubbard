#include "itensor/all.h"
#include "boson.h"
#include "TimeEvolve.hpp"
#include <time.h>


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

  //
  // Set the initial wavefunction matrix product state
  //
  auto state = InitState(sites);
  int p = Npart;
  for(int i = N; i >= 1; --i) {
      if (p >= 1) { state.set(i,"Occ1"); p -= 1; }
      else { state.set(i,"Emp"); }
  }

  auto psi = IQMPS(state);
  //Save initial state;
  auto psi0 = psi;


  // Regular time evolve
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

  auto psi1     = psi;
  auto args     = Args("Cutoff=",1E-8,"Maxm=",300);
  clock_t begin = clock();
  auto psi_t1   = TimeEvolve(psi1,ampolist,Cplx_i*tstep,args);
  clock_t end   = clock();
  std::cout << "Runtime for regular time evolution = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";
  printfln("Maximum MPS bond dimension after time evolution is %d",maxM(psi));
  Print(overlapC(psi_t1.back(),psi0));

          //Define the type "Gate" as a shorthand for BondGate<ITensor>
          using Gate = BondGate<IQTensor>;

          //Create a std::vector (dynamically sizeable array)
          //to hold the Trotter gates
          //
          auto gatelist = vector<Gate>();

            //Create the gates exp(-i*tstep/2*hterm)
            //and add them to gates
            for(int i = 1; i < N; ++i)
                {
                auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
                hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);

                // if (i%2) {
                //   hterm += 2.0*sites.op("N(N-1)",i)*sites.op("Id",i+1);
                //   hterm += 2.0*sites.op("Id",i)*sites.op("N(N-1)",i+1);
                // }

                auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
                gatelist.push_back(g);
                }
            //Create the gates exp(-i*tstep/2*hterm) in reverse
            //order (to get a second order Trotter breakup which
            //does a time step of "tstep") and add them to gates
            for(int i = N-1; i >= 1; --i)
                {
                auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
                hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);

                auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
                gatelist.push_back(g);
                }

          auto expU = vector<IQTensor>();
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
            expU.push_back(gate_);
          }

  //
  // ATTEMPT TO SPEED UP GATE CREATION!!
  //
  // using Gate = BondGate<IQTensor>;
  // std::vector< vector<Gate> > gatelist;
  // // First, load all operator tensors needed
  // std::vector<IQTensor> AAdaglistP;
  // std::vector<IQTensor> AdagAlistP;
  // std::vector<IQTensor> AAdaglistM;
  // std::vector<IQTensor> AdagAlistM;
  // for (int i = 1; i < N; ++i) {
  //   int j = N-i;
  //   AAdaglistP.emplace_back( sites.op("A",i)*sites.op("Adag",i+1) );
  //   AdagAlistP.emplace_back( sites.op("Adag",i)*sites.op("A",i+1) );
  //   AAdaglistM.emplace_back( sites.op("A",j)*sites.op("Adag",j+1) );
  //   AdagAlistM.emplace_back( sites.op("Adag",j)*sites.op("A",j+1) );
  // }
  // // Next, build gates for each TimeStepNum
  // for (size_t j = 0; j < ttotal/tstep; j++) {
  //   auto gates = vector<Gate>();
  //
  //   //Create the gates exp(-i*tstep/2*hterm)
  //   //and add them to gates
  //   for(int i = 1; i < N; ++i) {
  //     int k       = N-i;
  //
  //     auto htermP = (-1.0-tstep*j)*AAdaglistP[i-1];
  //     htermP     += (-1.0-tstep*j)*AdagAlistP[i-1];
  //     auto gP     = Gate(sites,i,i+1,Gate::tReal,tstep/2.,htermP);
  //
  //     auto htermM = (-1.0-tstep*j)*AAdaglistM[i-1];
  //     htermM     += (-1.0-tstep*j)*AdagAlistM[i-1];
  //     auto gM     = Gate(sites,k,k+1,Gate::tReal,tstep/2.,htermM);
  //
  //     gates.push_back(gP);
  //     gates.push_back(gM);
  //   }
  //   std::cout << "Built gate " << j << " out of " << ttotal/tstep << '\n';
  //   gatelist.push_back(gates);
  // }

  //
  //  ATTEMPT TO SPEED UP TIME EVOLVE BY USING MIXED MPO AND GATE DESCRIPTION
  //
            // using Gate = BondGate<IQTensor>;
            // using MPO  = MPOt<IQTensor>;
            //
            // std::vector<MPO> MPOlist;
            // std::vector<Gate> gatelist;
            //
            // // for (size_t j = 0; j < ttotal/tstep; j++) {
            // //   auto ampo = AutoMPO(sites);
            // //   for(int i = 1; i <= N; ++i) {
            // //     double U = (1.0+4.0*tstep*j);
            // //     // double U = 0;
            // //     ampo += U,"N(N-1)",i;
            // //   }
            // //   auto expH = toExpH<IQTensor>(ampo,tstep*Cplx_i);
            // //   MPOlist.push_back(expH);
            // // }
            // auto ampo = AutoMPO(sites);
            // for(int i = 1; i <= N; ++i) {
            //   // double U = (1.0+4.0*tstep*j);
            //   double U = 1;
            //   ampo += U,"N(N-1)",i;
            // }
            // auto expU = toExpH<IQTensor>(ampo,tstep*Cplx_i);
            //
            // for(int i = 1; i < N; ++i)
            //   {
            //   auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
            //   hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);
            //
            //   auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
            //   gatelist.push_back(g);
            //   }
            // for(int i = N-1; i >= 1; --i)
            //   {
            //   auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
            //   hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);
            //
            //   auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
            //   gatelist.push_back(g);
            //   }

  //Time evolve, overwriting psi when done
  begin = clock();
  // auto psi_t2 = TimeEvolveMixed(gatelist,MPOlist,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});
  auto psi_t2 = TimeEvolveMixed2(gatelist,expU,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});
  end = clock();
  std::cout << "Runtime for gate time evolution = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";
  printfln("Maximum MPS bond dimension after time evolution is %d",maxM(psi));

  //Print overlap of final state with initial state
  //(Will be complex so using overlapC which can return complex);
  Print(overlapC(psi,psi0));

  // Compare the evolution
  for (size_t i = 0; i < psi_t1.size(); i++) {
    Print(overlapC(psi_t1[i],psi_t2[i]));
  }


  return 0;
}
