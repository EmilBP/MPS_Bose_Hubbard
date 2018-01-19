#include "itensor/all.h"
#include "boson.h"
#include "TimeEvolve.hpp"
#include <time.h>


using namespace itensor;
using std::vector;

int main() {
  int N = 15;
  int Npart = 15;
  int HilbertDim = 5;
  Real tstep = 0.01; //time step (smaller is generally more accurate)
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


  // Create AutoMPO
  auto ampo = AutoMPO(sites);
  auto ampo2 = ampo;
  for(int i = 1; i < N; ++i) {
    ampo += -1.0,"A",i,"Adag",i+1;
    ampo += -1.0,"Adag",i,"A",i+1;
  }

  //Define the type "Gate" as a shorthand for BondGate<ITensor>
  using Gate = BondGate<IQTensor>;

  //Create a std::vector (dynamically sizeable array)
  //to hold the Trotter gates
  auto gates = vector<Gate>();

  //Create the gates exp(-i*tstep/2*hterm)
  //and add them to gates
  for(int i = 1; i < N; ++i)
      {
      auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);

      auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      gates.push_back(g);
      }
  //Create the gates exp(-i*tstep/2*hterm) in reverse
  //order (to get a second order Trotter breakup which
  //does a time step of "tstep") and add them to gates
  for(int i = N-1; i >= 1; --i)
      {
      auto hterm = -1.0*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -1.0*sites.op("Adag",i)*sites.op("A",i+1);

      auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      gates.push_back(g);
      }

  //Save initial state;
  auto psi0 = psi;


  // Regular time evolve
  auto psi1 = psi;
  auto args   = Args("Cutoff=",1E-8,"Maxm=",300);
  int Nsteps = ttotal/tstep;
  clock_t begin = clock();
  auto psi_t1 = TimeEvolve(psi1,ampo,Cplx_i*tstep,Nsteps,args);
  clock_t end = clock();
  std::cout << "Runtime for regular time evolution = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";
  printfln("Maximum MPS bond dimension after time evolution is %d",maxM(psi));
  Print(overlapC(psi_t1.back(),psi0));


  //Time evolve, overwriting psi when done
  begin = clock();
  auto psi_t2 = TimeEvolve(gates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});
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
