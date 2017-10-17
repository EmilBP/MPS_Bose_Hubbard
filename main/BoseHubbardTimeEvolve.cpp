  #include "itensor/all.h"
  #include "boson.h"
  #include "correlations.h"
  #include <vector>
  #include <fstream>

  using namespace itensor;
  using std::vector;

  int main()
    {
    int N = 50;
    int Npart = 50;

    auto sites = Boson(N);

    //
    // Create Super-Fluid (SF) Hamiltonian using AutoMPO
    //
    double J = -1.0;
    double U = 1e-2;
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

    auto H_SF = IQMPO(ampo);


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


    //
    // Set the parameters controlling the accuracy of the DMRG
    //
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,50,50,100,200;
    sweeps.cutoff() = 1E-9;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H_SF,sweeps,"Quiet");

    //
    // Save condensate fraction
    //
    vector<double> condensateFraction;

    auto rho = correlations::correlationTerm(sites,psi,"Adag","A");
    double fc = rho/Npart;
    condensateFraction.push_back(fc);


    //
    // Setup time evolution with Mott-Insulator (MI) Hamiltonian
    //
    J = -1e-4;
    U = 1;

    auto ampo2 = AutoMPO(sites);
    for(int i = 1; i < N; ++i) {
      ampo2 += J,"A",i,"Adag",i+1;
      ampo2 += J,"Adag",i,"A",i+1;
    }
    for (int i = 1; i <= N; ++i) {
      ampo2 += U/2.0,"N",i,"N-1",i;
      ampo2 += eps,"N",i;
    }

    auto H_MI = IQMPO(ampo2);
    auto tau = 1e-1;
    auto expH = toExpH<ITensor>(ampo2,tau*Cplx_i);


    //
    // Perform time evolution
    //
    vector<double> tvec = {0};
    auto args = Args("Cutoff=",1E-9,"Maxm=",200);
    auto ttotal = 8.0;
    auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));

    for(int n = 1; n <= nt; ++n){
        psi = exactApplyMPO(expH,psi,args);
        if ( (n*tau)%0.2 == 0) {
          rho = correlations::correlationTerm(sites,psi,"Adag","A");
          fc = rho/Npart;
          condensateFraction.push_back(fc);
          tvec.push_back(n*tau);
        }
    }

    std::fstream myfile;
    myfile.open("TimeEvolutionFractionData.txt",std::fstream::out);

    for (size_t i = 0; i < tvec.size(); i++) {
      myfile << tvec.at(i) << "\t";
      myfile << condensateFraction.at(i) << "\n";
    }
    myfile.close();


    return 0;
    }
