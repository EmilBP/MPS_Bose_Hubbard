#ifndef BOSEHUBBARDMPO_H
#define BOSEHUBBARDMPO_H

#include "autoMPOstrategy.h"
#include "boson.h"

using namespace itensor;

class BoseHubbardMPO : public autoMPOstrategy{
private:
  AutoMPO baseMPO;
  size_t N;

public:
  BoseHubbardMPO(SiteSet& sites, double J, double U, double eps);

  virtual AutoMPO updateMPO(double control);
  virtual AutoMPO derivative_control(double control);
};

#endif
