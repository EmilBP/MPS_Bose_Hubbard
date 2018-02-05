#ifndef HAMILTONIANBH_H
#define HAMILTONIANBH_H

#include "itensor/all.h"

using namespace itensor;

class HamiltonianBH {
private:
  SiteSet sites;
  double J;
  int N;

public:
  HamiltonianBH(const SiteSet& sites, const double J);

  IQMPO dHdU(const double& control_n, const double& tstep);
};

#endif
