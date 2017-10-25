#include "itensor/all.h"
#include "boson.h"
#include "correlations.h"
#include <vector>
#include <fstream>

using namespace itensor;
using std::vector;

int main(){

  // vector<int> NArray = {10,20,30,40,50};
  vector<int> NArray = {50};
  int Nnumber = NArray.size();

  //
  // Set sweep settinngs
  //
  auto sweeps = Sweeps(5); // set to min 3.
  sweeps.maxm() = 10,20,50,100,100,100,200;
  sweeps.cutoff() = 1E-9;
  sweeps.niter() = 2;
  sweeps.noise() = 1E-7,1E-8,0.0;
  println(sweeps);

  //
  // Set values for U
  //
  int Nsteps = 20; //set to 100.
  double Umin = 0, Umax = 0.1, step = Umax/Nsteps;
  vector<double> array;
  while(Umin <= Umax) {
      array.push_back(Umin);
      Umin += step;
  }

  double output[Nsteps][Nnumber+1];

  for (int j = 0; j < Nnumber; j++) {
    int N = NArray.at(j);
    int Npart = N;
    int HilbertDim = 10;
    auto sites = Boson(N,HilbertDim);
    //
    // Set the initial wavefunction matrix product state
    //
    auto state = InitState(sites);
    int p = Npart;
    for(int i = N; i >= 1; --i)
        {
        if (p >= 1) {
          state.set(i,"Occ1");
          p -= 1;
        }
        else
            {
            state.set(i,"Emp");
            }
        }

    auto psi = IQMPS(state);


    for (size_t Ui = 0; Ui < Nsteps; ++Ui) {

      double J = -1.0;
      double U = array.at(Ui);
      double eps = 0;

      //
      // build Hamiltonian MPO
      //
      auto ampo = AutoMPO(sites);
      for(int i = 1; i < N; ++i) {
        ampo += J,"A",i,"Adag",i+1;
        ampo += J,"Adag",i,"A",i+1;
      }
      for (int i = 1; i <= N; ++i) {
        ampo += U/2.0,"N",i,"N-1",i;
        ampo += eps,"N",i;
      }
      auto H = IQMPO(ampo);

      //
      //  Optimize with DMRG and get condensate fraction
      //
      auto energy = dmrg(psi,H,sweeps,"Quiet");
      auto lambda1 = correlations::correlationTerm(sites,psi,"Adag","A");
      double fc = lambda1/Npart;
      output[Ui][j+1] = fc;
      output[Ui][0] = U;
    }
  }

  std::fstream myfile;
  myfile.open("condensateData2.txt",std::fstream::out);

  for (size_t i = 0; i < Nsteps; i++) {
    for (size_t j = 0; j < Nnumber+1; j++) {
      myfile << output[i][j] << "\t";
    }
    myfile << "\n";
  }
  myfile.close();

  return 0;
}
