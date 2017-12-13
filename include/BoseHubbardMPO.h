#ifndef BOSEHUBBARDMPO_H
#define BOSEHUBBARDMPO_H

#include "autoMPOstrategy.h"
#include "boson.h"
#include <math.h>

using namespace itensor;

class BoseHubbardMPO : public autoMPOstrategy{
private:
  AutoMPO baseMPO;
  size_t N;
  double trapStr;
  int L;

public:
  BoseHubbardMPO(SiteSet& sites, double trapStr = 0);
  double getJ(double V0);
  double getU(double V0);
  double getdJdV(double V0);
  double getdUdV(double V0);

  virtual AutoMPO updateMPO(double control);
  virtual AutoMPO dHdu(double control);
};

#endif
