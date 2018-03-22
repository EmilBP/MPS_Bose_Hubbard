#ifndef OPTIMALCONTROL_H
#define OPTIMALCONTROL_H

#include "itensor/all.h"
#include "ControlBasis.hpp"

#include <vector>
#include <iterator>
#include <algorithm>



using namespace itensor;
using vec = std::vector<double>;
using vecpair = std::pair<double, vec>;

template<class TimeStepper, class Hamiltonian>
class OptimalControl{
private:
  TimeStepper& timeStepper;
  Hamiltonian& hamil;
  double gamma, tstep;
  IQMPS psi_target, psi_init;

  std::vector<IQMPS> psi_t;
  std::vector<IQMPS> chi_t;

  double getFidelity(const vec& control);
  double getRegularisation(const vec& control);
  vecpair getRegPlusRegGrad(const vec& control);
  vecpair getFidelityPlusFidelityGrad(const vec& control);
  void calcPsi(const vec& control);
  void calcChi(const vec& control);

public:
  OptimalControl(IQMPS& psi_target, IQMPS& psi_init, TimeStepper& timeStepper, Hamiltonian& hamil, double gamma);

  double getCost(const vec& control);
  vecpair getAnalyticGradient(const vec& control);
  vecpair getNumericGradient(const vec& control);

  double getCost(const ControlBasis& bControl);
  vecpair getAnalyticGradient(const ControlBasis& bControl);
  vecpair getNumericGradient(const ControlBasis& bControl);
  vec getFidelityForAllT(const vec& control);
  vec getFidelityForAllT(const ControlBasis& bControl);

  vecpair checkCostPlusFidelity(const vec& control);
};

#endif
