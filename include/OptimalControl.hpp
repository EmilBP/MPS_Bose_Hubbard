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

class OptimalControl{
private:
  double gamma;
  IQMPS psi_target, psi_init;
  autoMPOstrategy& MPOstrat;

  std::vector<IQMPS> psi_t;
  std::vector<IQMPS> chi_t;

  double getCost(IQMPS& psi, std::vector<double>& control);
  std::vector<AutoMPO> updateHamiltonian(std::vector<double>& control);
  std::vector<double> updateControl(std::vector<double>& gradient);

public:
  OptimalControl(IQMPS& psi_target, IQMPS& psi_init, autoMPOstrategy& MPOstrat, double gamma);
  std::vector<double> Optimize(std::vector<double>& control_init, double dt, size_t maxeval, const Args& args);

  std::vector<double> getAnalyticGradient(std::vector<double>& control,double dt);
  std::vector<double> getAnalyticGradientTest(std::vector<double>& control,double dt,const Args& args);
  std::vector<double> getNumericGradient(std::vector<double>& control,double epsilon, double dt, const Args& args);
};

#endif
