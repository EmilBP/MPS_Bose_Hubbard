#include "OptimalControl.hpp"

OptimalControl::CostFunction(IQMPS& psi, std::vector<double>& control){
  // mssing gamma term
  return 0.5*(1- overlap(psi_target,psi));
}

std::vector<AutoMPO> OptimalControl::updateMPO(AutoMPO (*func)(double),std::vector<double>& control){
  std::vector<AutoMPO> v;
  v.reserve(size(control));
  for (auto &u : control) {
    v.emplace_back(baseMPO + (*func)(u));
  }
  return v;
}

std::vector<double> OptimalControl::OptimizeControl(std::vector<double>& control_init, size_t maxeval){
  IQMPS chi_T;
  std::vector<double> gradient;
  std::vector<AutoMPO> Ampos;


  for (size_t i = 0; i < maxeval; i++) {
    Ampos = updateMPO(func, control_init);

    psi_t = TimeEvolve(psi_init, ampo, tau, args);
    chi_T = -Cplx_i * psi_target*overlap(psi_target,psi_t.back());
    chi_t = TimeEvolveBack(chi_t, ampo, tau, args);

    gradient = calculateGradient();
    control_init = updateControl(gradient);
  }
}
