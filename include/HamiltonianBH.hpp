#ifndef HAMILTONIANBH_H
#define HAMILTONIANBH_H

#include "itensor/all.h"
#include <complex>

using namespace itensor;

class HamiltonianBH {
private:
  SiteSet sites;
  double J, tstep;
  int N;
  size_t expansionOrder;
  std::vector<std::complex <double> > prefactors;
  AutoMPO ampo;
  IQMPO dHdUconst;


public:
  HamiltonianBH(const SiteSet& sites, double J, double tstep, size_t expansionOrder = 0);

  void setTstep(const double tstep_);
  void setExpansionOrder(const size_t expansionOrder_);
  IQMPO dHdU(const double& control_n);
};

#endif
