#ifndef OPTIMALCONTROL_H
#define OPTIMALCONTROL_H

#include "itensor/all.h"
#include "TimeEvolve.hpp"
#include "autoMPOstrategy.h"
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>


namespace itensor{

class OptimalControl{
private:
  double gamma;
  IQMPS psi_target, psi_init;
  autoMPOstrategy& MPOstrat;

  std::vector<IQMPS> psi_t;
  std::vector<IQMPS> chi_t;

  double CostFunction(IQMPS& psi, std::vector<double>& control);
  std::vector<AutoMPO> updateMPOlist(std::vector<double>& control);
  std::vector<double> calculateGradient(std::vector<double>& control,double dt);
  std::vector<double> updateControl(std::vector<double>& gradient);

public:
  OptimalControl(IQMPS& psi_target, IQMPS& psi_init, autoMPOstrategy& MPOstrat, double gamma);
  std::vector<double> OptimizeControl(std::vector<double>& control_init, double dt, size_t maxeval, const Args& args);

};

} // end namespace
#endif
