#ifndef OPTIMALCONTROL_H
#define OPTIMALCONTROL_H

#include "itensor/all.h"
#include <vector>
#include <iterator>
#include <algorithm>



using namespace itensor;
using vec = const std::vector<double>;

template<class TimeStepper>
class OptimalControl{
private:
  double gamma, tstep;
  IQMPS psi_target, psi_init;
  TimeStepper timeStepper;
  MPOt<IQTensor> dHdU;

  std::vector<IQMPS> psi_t;
  std::vector<IQMPS> chi_t;

  double getFidelity(vec& control);
  double getRegularisation(vec& control);
  std::vector<double> getRegGrad(vec& control);
  std::vector<double> getFidelityGrad(vec& control);
  void calcPsi(vec& control);
  void calcChi(vec& control);


public:
  OptimalControl(IQMPS& psi_target, IQMPS& psi_init, TimeStepper& timeStepper, MPOt<IQTensor>& dHdU, double gamma, double tstep);

  double getCost(vec& control);
  std::vector<double> getAnalyticGradient(vec& control);
  std::vector<double> getNumericGradient(vec& control);
};

#endif
