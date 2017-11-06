#include "BoseHubbardMPO.h"

BoseHubbardMPO::BoseHubbardMPO(AutoMPO baseMPO, size_t N)
  : baseMPO(baseMPO), N(N) {

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
