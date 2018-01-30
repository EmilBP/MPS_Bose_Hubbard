#ifndef OPTIMALCONTROL_H
#define OPTIMALCONTROL_H

#include "itensor/all.h"
#include "TimeEvolve.hpp"
#include "autoMPOstrategy.h"
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>


using namespace itensor;

template<class TimeStepper>
class OptimalControl{
private:
  double gamma, tstep;
  IQMPS psi_target, psi_init;
  TimeStepper timeStepper;

  std::vector<IQMPS> psi_t;
  std::vector<IQMPS> chi_t;

  // void initGates();
  double getFidelity(const std::vector<double>& control);
  double getRegularisation(const std::vector<double>& control);
  std::vector<double> getRegGrad(const std::vector<double>& control);
  std::vector<double> getFidelityGrad(const std::vector<double>& control);


public:
  OptimalControl(IQMPS& psi_target, IQMPS& psi_init, TimeStepper& timeStepper, double gamma, double tstep);

  double getCost(const std::vector<double>& control);
  std::vector<double> getAnalyticGradient(const std::vector<double>& control);
  std::vector<double> getNumericGradient(const std::vector<double>& control);
  // std::vector<double> getAnalyticGradientTest(std::vector<double>& control,const Args& args);
};

#endif
