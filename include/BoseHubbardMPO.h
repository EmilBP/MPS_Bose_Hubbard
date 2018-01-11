#ifndef BOSEHUBBARDMPO_H
#define BOSEHUBBARDMPO_H

#include "autoMPOstrategy.h"
#include "boson.h"
#include <math.h>
#include <string>

using namespace itensor;

typedef std::vector<std::vector<double>> matrix;

class BoseHubbardMPO : public autoMPOstrategy{
private:
  AutoMPO baseMPO;
  size_t N;
  double overlayTrapDepth;
  int L;
  matrix UJdata;
  bool hasData;

public:
  BoseHubbardMPO(SiteSet& sites, double overlayTrapDepth = 0);
  void interpolateData(double V0, double& U, double& J, double& dUdV, double& dJdV);
  void loadUJdata(std::string& filename);

  virtual AutoMPO updateMPO(double control);
  virtual AutoMPO dHdu(double control);
};

#endif
