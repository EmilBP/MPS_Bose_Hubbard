#ifndef TIMESTEPPERTEBD_H
#define TIMESTEPPERTEBD_H

#include "itensor/all.h"

using namespace itensor;

class TimeStepperTEBD {
private:
  double tstep, J;
  std::vector< BondGate<IQTensor> > JGates;
  std::vector<IQTensor> UGates;
  SiteSet sites;
  Args args;

  void initJGates(const double J);
  void initUGates(const double U);

public:
  TimeStepperTEBD(const SiteSet& sites, const double J, const double tstep, const Args& args);
  void setTstep(const double tstep_);
  void step(IQMPS& psi, const double from, const double to);
};

#endif
