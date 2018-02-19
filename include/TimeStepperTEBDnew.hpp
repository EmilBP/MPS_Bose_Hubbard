#ifndef TIMESTEPPERTEBDNEW_H
#define TIMESTEPPERTEBDNEW_H

#include "itensor/all.h"

using namespace itensor;
using GateList = std::vector< BondGate<IQTensor> >;

class TimeStepperTEBDnew {
private:
  double tstep, J;
  GateList JGates_tforwards;
  GateList JGates_tbackwards;
  std::vector<IQTensor> UGates;
  SiteSet sites;
  Args args;

  void initJGates(const double J);
  void initUGates(const double U);
  void doStep(IQMPS& psi, const GateList& JGates);

public:
  TimeStepperTEBDnew(const SiteSet& sites, const double J, const double tstep, const Args& args);
  void setTstep(const double tstep_);
  void step(IQMPS& psi, const double from, const double to, bool propagateForward = true);
  double getTstep();
};

#endif
