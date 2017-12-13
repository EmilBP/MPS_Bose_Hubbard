#include "BoseHubbardMPO.h"

BoseHubbardMPO::BoseHubbardMPO(SiteSet& sites, double trapStr)
  : baseMPO(AutoMPO(sites)), N(sites.N()), trapStr(trapStr), L(sites.N()) {

}

double BoseHubbardMPO::getJ(double V0){
  // in units of E_rec
  return 4.0/sqrt(M_PI)*pow(V0,0.75)*exp(-2.0*sqrt(V0));
}

double BoseHubbardMPO::getU(double V0){
  // in units of E_rec
  double ka = 1;
  return sqrt(8.0/M_PI)*ka*pow(V0,0.75);
}

double BoseHubbardMPO::getdJdV(double V0){
  return exp(-2.0*sqrt(V0)) * (1.69257 - 2.25676*sqrt(V0))*pow(V0,-0.25);
}

double BoseHubbardMPO::getdUdV(double V0){
  double ka = 1;
  return (1.19683*ka)*pow(V0,-0.25);
}

AutoMPO BoseHubbardMPO::updateMPO(double control){
  AutoMPO ampo  = baseMPO;
  double J      = getJ(control);
  double U      = getU(control);

  for(int i = 1; i < N; ++i) {
    ampo += -J,"A",i,"Adag",i+1;
    ampo += -J,"Adag",i,"A",i+1;
  }
  for (int i = 1; i <= N; ++i) {
    ampo += U/2.0,"N(N-1)",i;
    ampo += 0.5*trapStr*(i-0.5*(L-1))*(i-0.5*(L-1)),"N",i;
  }

  return ampo;
}

AutoMPO BoseHubbardMPO::dHdu(double control){
  AutoMPO ampo = baseMPO;
  double dJdV  = getdJdV(control);
  double dUdV  = getdUdV(control);

  for(int i = 1; i < N; ++i) {
    ampo += -dJdV,"A",i,"Adag",i+1;
    ampo += -dJdV,"Adag",i,"A",i+1;
  }
  for (int i = 1; i <= N; ++i) {
    ampo += dUdV/2.0,"N(N-1)",i;
  }

  return ampo;
}
