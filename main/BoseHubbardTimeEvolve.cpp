  #include "itensor/all.h"
  #include "boson.h"
  #include "correlations.h"
  #include "TimeEvolve.hpp"
  #include <vector>
  #include <fstream>

  using namespace itensor;
  using std::vector;

  int main()
    {
    int N = 5;
    int Npart = 5;
    int HilbertDim = 5;

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

    //
    // Create Super-Fluid (SF) Hamiltonian using AutoMPO
    //
    double J = -1.0;
    double U = 0;
    double eps = 0;

    auto ampo = AutoMPO(sites);
    auto ampo2 = ampo;
    for(int i = 1; i < N; ++i) {
      ampo += J,"A",i,"Adag",i+1;
      ampo += J,"Adag",i,"A",i+1;
    }
    for (int i = 1; i <= N; ++i) {
      ampo += U/2.0,"N(N-1)",i;
      ampo += eps,"N",i;
    }

    auto H_SF = IQMPO(ampo);


    //
    // Set the parameters controlling the accuracy of the DMRG
    //
    auto sweeps = Sweeps(10);
    sweeps.maxm() = 10,20,50,50,100,200;
    sweeps.cutoff() = 1E-9;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H_SF,sweeps,"Quiet");
    printf("Energy w.r.t. superfluid Hamiltonian: %f\n",energy);
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
    auto tau = 1e-2;

    J = 0;
    U = 1;

    for(int i = 1; i < N; ++i) {
      ampo2 += J,"A",i,"Adag",i+1;
      ampo2 += J,"Adag",i,"A",i+1;
    }
    for (int i = 1; i <= N; ++i) {
      ampo2 += U/2.0,"N(N-1)",i;
    }

    auto H_MI = IQMPO(ampo2);
//    auto expH = toExpH<IQTensor>(ampo2,tau*Cplx_i);
    auto expiHdt_1 = toExpH<IQTensor>(ampo2,tau*0.5*(1.0+Cplx_i)*Cplx_i);
    auto expiHdt_2 = toExpH<IQTensor>(ampo2,tau*0.5*(1.0-Cplx_i)*Cplx_i);

    //
    // Perform time evolution
    //
    vector<double> tvec = {0};
    auto args = Args("Cutoff=",1E-9,"Maxm=",50);
    auto ttotal = 3;
    auto nt = int(ttotal/tau);


    for(int n = 1; n <= nt; ++n){
//        psi = exactApplyMPO(expH,psi,args);
//        normalize(psi);
      fitApplyMPO(psi,expiHdt_1,psi,args);
      fitApplyMPO(psi,expiHdt_2,psi,args);
	if (n % 10) {
            tvec.push_back(n*tau);
            auto rho = correlations::correlationTerm(sites,psi,"Adag","A");
            auto fc = rho/Npart;
            condensateFraction.push_back(fc);
	}
    }

    std::cout << condensateFraction.back() << std::endl;



    return 0;
    }
