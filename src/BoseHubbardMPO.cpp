#include "BoseHubbardMPO.h"

BoseHubbardMPO::BoseHubbardMPO(SiteSet& sites, double J, double U, double eps)
  : baseMPO(AutoMPO(sites)), N(sites.N()) {

  for(int i = 1; i < N; ++i) {
    baseMPO += J,"A",i,"Adag",i+1;
    baseMPO += J,"Adag",i,"A",i+1;
  }
  for (int i = 1; i <= N; ++i) {
    baseMPO += U/2.0,"N(N-1)",i;
    baseMPO += eps,"N",i;
  }
}

AutoMPO BoseHubbardMPO::updateMPO(double control){
  AutoMPO ampo = baseMPO;
  for (int i = 1; i <= N; ++i) {
    ampo += control/2.0,"N(N-1)",i;
  }
  return ampo;
}

AutoMPO BoseHubbardMPO::derivative_control(double control){
  AutoMPO ampo = baseMPO;
  for (int i = 1; i <= N; ++i) {
    ampo += 0.5,"N(N-1)",i;
  }
  return ampo;
}
