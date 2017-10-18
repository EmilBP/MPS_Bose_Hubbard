#include "itensor/all.h"
#include "boson.h"
#include "correlations.h"

using namespace itensor;

int main()
    {
    int N = 50;
    int Npart = 50;
    int HilbertDim = 20;

    auto sites = Boson(N,HilbertDim);

    //
    // Create the Hamiltonian using AutoMPO
    //
    double J = -1.0;
    double U = 0;
    double eps = 0;

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
    // Set the initial wavefunction matrix product state
    //
    auto state = InitState(sites);
    int p = Npart;
    for(int i = N; i >= 1; --i)
        {
        if (p >= 1) {
          println("Singly occupying site ",i);
          state.set(i,"Occ1");
          p -= 1;
        }
        else
            {
            state.set(i,"Emp");
            }
        }

    auto psi = IQMPS(state);

    printfln("Initial energy = %.5f", overlap(psi,H,psi) );

    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).

    auto sweeps = Sweeps(10);
    sweeps.maxm() = 20,30,50,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-6,1E-7,1E-8,1E-9,0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("Wavefunc norm = %.10f",norm(psi));

    auto test = correlations::correlationTerm(sites,psi,"Adag","A");
    double fc = test/Npart;
    printfln("\nCorrelation = %.10f",fc);

    return 0;
    }
