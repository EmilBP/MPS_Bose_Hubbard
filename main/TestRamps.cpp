#include "itensor/all.h"
#include "boson.h"
#include "OptimalControl.hpp"
#include "correlations.hpp"
#include "InitializeState.hpp"
#include "BoseHubbardMPO.h"
#include <math.h>

using namespace itensor;


double getRamp(double t, double T){
  return 3+11/T*t;
}

int main(){
  int N       = 32;
  int Npart   = 16;
  int locDim  = 8;
  auto sites  = Boson(N,locDim);

  //
  //  setup trap
  //
  double m        = 1;
  double omega    = 1;
  double a_lat    = 1;
  double trapStr  = m*omega*omega*a_lat*a_lat;
  auto BHMPO      = BoseHubbardMPO(sites,trapStr);

  //
  //  setup initial state
  //
  auto ampo_init  = BHMPO.updateMPO(3); // initial trap depth
  auto H_init     = IQMPO(ampo_init);
  auto state      = InitState(sites);
  int p           = Npart;
  for(int i = 3*N/4; i >= 1; --i){
      if (p >= 1) {
        println("Singly occupying site ",i);
        state.set(i,"Occ1");
        p -= 1;
      }
      else { state.set(i,"Emp"); }
  }
  auto psi        = IQMPS(state);

  auto sweeps     = Sweeps(5);
  sweeps.maxm()   = 20,30,50,100,200;
  sweeps.cutoff() = 1E-8;
  sweeps.niter()  = 2;
  sweeps.noise()  = 1E-7,1E-8,0;
  auto energy     = dmrg(psi,H_init,sweeps,"Quiet");

  //
  //  setup ramp
  //

  auto args   = Args("Cutoff=",1E-8,"Maxm=",200);
  double eps  = 1e-6;
  double T    = 11.75 * 1e-3;
  double dt   = T/1e2;
  double temp = 0;

  std::vector<double> times;
  std::vector<double> ramp;
  std::vector<AutoMPO> mpos;

  times.emplace_back(temp);
  ramp.emplace_back(getRamp(temp,T));
  while (T-temp > 1e-8) {
    temp += dt;
    times.emplace_back(temp);
    ramp.emplace_back(getRamp(temp,T));
  }

  for (auto& u : ramp){
    std::cout << "V0: " << u << " J: " << BHMPO.getJ(u) << ". U: " << BHMPO.getU(u) << std::endl;
    mpos.emplace_back(BHMPO.updateMPO(u));
  }

  auto psi_t = TimeEvolve(psi,mpos,Cplx_i*dt,args);

  //
  //  calculate figure of merit from Frank/Bloch article
  //
  std::vector<Cplx> variances;
  for (auto& psi_i : psi_t){
    Cplx vartemp = 0;
    for (size_t i = N/2-3; i < N/2+4; i++) {
      auto var0 = expectationValue(sites,psi,"NN",i) - pow(expectationValue(sites,psi,"N",i),2);
      vartemp  += (expectationValue(sites,psi_i,"NN",i) - pow(expectationValue(sites,psi_i,"N",i),2))/var0;
    }
    variances.emplace_back(vartemp/8.0);
  }

  //
  // print
  //
  std::cout << "Variances:" << std::endl;
  for (auto& v : variances){
    std::cout << "\t" << v << "\n";
  }

  return 0;
}
