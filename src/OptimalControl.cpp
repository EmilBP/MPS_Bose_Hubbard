#include "OptimalControl.hpp"

OptimalControl::OptimalControl(IQMPS& psi_target, IQMPS& psi_init, autoMPOstrategy& MPOstrat, double gamma)
  : psi_target(psi_target), psi_init(psi_init), MPOstrat(MPOstrat), gamma(gamma){

}

double OptimalControl::getCost(IQMPS& psi, std::vector<double>& control){
  // mssing gamma term
  Cplx proj = overlap(psi_target,psi);
  double normsquare = Real(proj*conj(proj));
  return 0.5*(1-normsquare);
}

std::vector<AutoMPO> OptimalControl::updateHamiltonian(std::vector<double>& control){
  std::vector<AutoMPO> v;
  v.reserve(control.size());
  for (auto &u : control) {
    v.emplace_back(MPOstrat.updateMPO(u));
  }
  return v;
}

std::vector<double> OptimalControl::getAnalyticGradient(std::vector<double>& control, double dt){
  // missing gamma term
  std::vector<double> g;
  g.reserve(control.size());

  auto i1 = begin(chi_t), i2 = begin(psi_t), e = end(chi_t);
  auto i3 = begin(control);

  while (i1 != e) {
    g.emplace_back(dt*Real(overlap(*i1,IQMPO(MPOstrat.dHdu(*i3)),*i2)));
    ++i1;
    ++i2;
    ++i3;
  }
  return g;
}

std::vector<double> OptimalControl::getNumericGradient(
      std::vector<double>& control, double epsilon,
      double dt, const Args& args)
{
  double Jp, Jm;
  std::vector<double> g;
  std::vector<IQMPS> psi_temp;
  std::vector<AutoMPO> tempMPO;
  g.reserve(control.size());

  for (auto& ui : control){
    ui        += epsilon;
    tempMPO    = updateHamiltonian(control);
    psi_temp   = TimeEvolve(psi_init,tempMPO,dt*Cplx_i,args);
    Jp         = getCost(psi_temp.back(),control);

    ui        -= 2*epsilon;
    tempMPO    = updateHamiltonian(control);
    psi_temp  = TimeEvolve(psi_init,tempMPO,dt*Cplx_i,args);
    Jm         = getCost(psi_temp.back(),control);

    ui += epsilon;
    g.emplace_back((Jp-Jm)/(2.0*epsilon));
  }
  return g;
}

std::vector<double> OptimalControl::updateControl(std::vector<double>& gradient){
  // PLACEHOLDER
  return gradient;
}

std::vector<double> OptimalControl::Optimize(
      std::vector<double>& control_init,
      double dt,
      size_t maxeval,
      const Args& args)
{
  IQMPS chi_T;
  std::vector<double> gradient;
  std::vector<AutoMPO> AmpoList;
  std::vector<double> costs;
  costs.reserve(maxeval+1);

  for (size_t i = 0; i < maxeval; i++) {
    AmpoList = updateHamiltonian(control_init);

    psi_t = TimeEvolve(psi_init, AmpoList, dt*Cplx_i, args);
    chi_T = -Cplx_i * psi_target*overlap(psi_target,psi_t.back());
    chi_t = TimeEvolveBack(chi_T, AmpoList, dt*Cplx_i, args);

    gradient = getAnalyticGradient(control_init,dt);
    control_init = updateControl(gradient);
    costs.emplace_back(getCost(psi_t.back(),control_init));
  }
  //
  //  get final cost
  //
  AmpoList = updateHamiltonian(control_init);
  psi_t = TimeEvolve(psi_init, AmpoList, dt*Cplx_i, args);
  costs.emplace_back(getCost(psi_t.back(),control_init));

  return costs;
}
