#ifndef TIMESTEPPERTEBDFAST_H
#define TIMESTEPPERTEBDFAST_H

#include "itensor/all.h"

using namespace itensor;
using GateList = std::vector< BondGate<IQTensor> >;

class TimeStepperTEBDfast {
private:
  double tstep, J;
  GateList JGates_forwards;
  GateList JGates_backwards;
  std::vector<IQTensor> UGates;
  SiteSet sites;
  Args args;

  void initJGates(const double J);
  void initUGates(const double U);
  void doStep(IQMPS& psi, const GateList& JGates, const std::vector<IQTensor> UGates);

public:
  TimeStepperTEBDfast(const SiteSet& sites, const double J, const double tstep, const Args& args);
  void setTstep(const double tstep_);
  void step(IQMPS& psi, const double from, const double to, bool propagateForward = true);
  double getTstep();
};

#endif
