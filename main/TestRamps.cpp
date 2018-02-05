#include "itensor/all.h"
#include "boson.h"
#include "OptimalControl.hpp"
#include "TimeEvolve.hpp"
#include "correlations.hpp"
#include "InitializeState.hpp"
#include <math.h>
#include <string>
#include <fstream>


using namespace itensor;


double getRamp(double t, double T){
  return 3+11/T*t;
}

int main(){
  // int N       = 32;
  // int Npart   = 16;
  // int locDim  = 8;
  // auto sites  = Boson(N,locDim);
  //
  // //
  // //  setup trap
  // //
  // double trapStr  = 0.00242; // 0.5*m*omega^2 a_lat² in units of E_rec
  // auto BHMPO      = BoseHubbardMPO(sites,trapStr);
  // std::string filename = "UJparams";
  // BHMPO.loadUJdata(filename);
  //
  // //
  // //  setup initial state
  // //
  // auto ampo_init  = BHMPO.updateMPO(3); // initial trap depth
  // auto H_init     = IQMPO(ampo_init);
  // auto state      = InitState(sites);
  // int p           = Npart;
  // for(int i = N; i >= 1; --i){
  //     if (i <= 3*N/4 && p >= 1) {
  //       state.set(i,"Occ1");
  //       p -= 1;
  //     }
  //     else {
  //       state.set(i,"Emp");
  //     }
  // }
  // // for(int i = N; i >= 1; --i){
  // //     if ( p >= 1) {
  // //       state.set(i,"Occ1");
  // //       p -= 1;
  // //     }
  // //     else {
  // //       state.set(i,"Emp");
  // //     }
  // // }
  // auto psi        = IQMPS(state);
  //
  // auto sweeps     = Sweeps(5);
  // sweeps.maxm()   = 20,30,50,100,200,250,300;
  // sweeps.cutoff() = 1E-9;
  // sweeps.niter()  = 2;
  // sweeps.noise()  = 1E-7,1E-8,0;
  // auto energy     = dmrg(psi,H_init,sweeps,"Quiet");
  //
  // auto psi_0 = psi;
  //
  // for (size_t i = 1; i <= N; i++) {
  //   std::cout << expectationValue(sites,psi_0,"N",i) << '\n';
  // }
  // //
  // //  setup ramp
  // //
  //
  // auto args   = Args("Cutoff=",1E-9,"Maxm=",300);
  // double eps  = 1e-6;
  // double T    = 11.75*1e-3 *12741.13; // t/hbar units of (E_rec)^-1
  // double dt   = T/1e2;
  // double temp = 0;
  //
  // std::vector<double> times;
  // std::vector<double> ramp;
  // std::vector<AutoMPO> mpos;
  //
  // times.emplace_back(temp);
  // ramp.emplace_back(getRamp(temp,T));
  // while (T-temp > 1e-8) {
  //   temp += dt;
  //   times.emplace_back(temp);
  //   ramp.emplace_back(getRamp(temp,T));
  // }
  //
  // for (auto& u : ramp){
  //   mpos.emplace_back(BHMPO.updateMPO(u));
  // }
  //
  // auto psi_t = TimeEvolve(psi,mpos,Cplx_i*dt,args);
  //
  // //
  // //  calculate figure of merit from Frank/Bloch article
  // //
  // std::vector<Cplx> variances;
  // for (auto& psi_i : psi_t){
  //   Cplx vartemp = 0;
  //   for (size_t i = N/2-4; i < N/2+4; i++) {
  //     double var0 = real(expectationValue(sites,psi_0,"NN",i) - pow(expectationValue(sites,psi_0,"N",i),2));
  //     vartemp  += (expectationValue(sites,psi_i,"NN",i) - pow(expectationValue(sites,psi_i,"N",i),2))/var0;
  //   }
  //   variances.emplace_back(vartemp/8.0);
  // }
  //
  // //
  // // print ans save data
  // //
  // // std::ofstream outFile("RampTestData.txt");
  //
  // // for (size_t i = 0; i < ramp.size(); i++) {
  // //   double U,J,temp,temp2;
  // //   BHMPO.interpolateData(ramp[i],U,J,temp,temp2);
  // //   std::cout << "t = " << times[i]*1e3 /12741.13 << '\t';
  // //   std::cout << "V0 = " << ramp[i] << '\t';
  // //   std::cout << "U = " << U << '\t';
  // //   std::cout << "J = " << J << '\t';
  // //   std::cout << "U/J = " << U/J << '\t';
  // //   std::cout << "F = " << real(variances[i]) << '\n';
  // //
  // //   outFile << times[i]*1e3 /12741.13 << '\t';
  // //   outFile << ramp[i] << '\t';
  // //   outFile << U << '\t';
  // //   outFile << J << '\t';
  // //   outFile << U/J << '\t';
  // //   outFile << real(variances[i]) << '\n';
  // // }
  //
  // for (size_t i = 1; i <= N; i++) {
  //   std::cout << expectationValue(sites,psi_t.back(),"N",i) << '\n';
  // }

  return 0;
}
