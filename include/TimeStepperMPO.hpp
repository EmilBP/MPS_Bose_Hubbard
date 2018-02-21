#ifndef TIMESTEPPERMPO_H
#define TIMESTEPPERMPO_H

#include "itensor/all.h"

using namespace itensor;

using splitstep = std::pair<IQMPO,IQMPO>;

class TimeStepperMPO {
private:
  double tstep, J;
  AutoMPO Jampo, ampo;
  SiteSet sites;
  Args args;

  void initJampo(const double J);
  void initUampo(const double U);
  void doStep(IQMPS& psi, const splitstep& tOp);

public:
  TimeStepperMPO(const SiteSet& sites, const double J, const double tstep, const Args& args);
  void setTstep(const double tstep_);
  void step(IQMPS& psi, const double from, const double to, bool propagateForward = true);
  double getTstep();
};

#endif
