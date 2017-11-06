#ifndef BOSEHUBBARDMPO_H
#define BOSEHUBBARDMPO_H

#include "autoMPOstrategy.h"
#include "boson.h"

class BoseHubbardMPO : public autoMPOstrategy{
private:
  AutoMPO baseMPO;
  size_t N;

public:
  BoseHubbardMPO(AutoMPO baseMPO, size_t N);

  virtual AutoMPO updateMPO(double control);
  virtual AutoMPO derivative_control(double control);
};

#endif
