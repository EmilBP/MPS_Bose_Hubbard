#ifndef OPTIMALCONTROL_H
#define OPTIMALCONTROL_H

#include "itensor/all.h"
#include "TimeEvolve.hpp"
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>


namespace itensor{

class OptimalControl{
private:
  double gamma;
  IQMPS psi_target;
  AutoMPO baseMPO;

  std::vector<IQMPS> psi_t;
  std::vector<IQMPS> chi_t;

  double CostFunction(IQMPS& psi, std::vector<double>& control);
  std::vector<AutoMPO> updateMPO(AutoMPO (*func)(double),std::vector<double>& control);

public:

  OptimalControl(IQMPS& psi_target, AutoMPO& baseMPO);
  std::vector<double> OptimizeControl(std::vector<double>& control_init, size_t maxeval);

};

} // end namespace
#endif
