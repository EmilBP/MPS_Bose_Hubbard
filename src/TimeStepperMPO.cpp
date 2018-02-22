#include "TimeStepperMPO.hpp"

TimeStepperMPO::TimeStepperMPO(const SiteSet& sites, const double J, const double tstep, const Args& args)
  : J(J), sites(sites), args(args), Jampo(AutoMPO(sites)), ampo(AutoMPO(sites)) {
  setTstep(tstep);
}

void TimeStepperMPO::initJampo(const double J){
  Jampo = AutoMPO(sites);

  for(int i = 1; i < sites.N(); ++i) {
    Jampo += -J,"Adag",i,"A",i+1;
    Jampo += -J,"A",i,"Adag",i+1;
  }
}


void TimeStepperMPO::setTstep(const double tstep_){
  tstep = tstep_;
  initJampo(J);
}

double TimeStepperMPO::getTstep(){
  return tstep;
}

void TimeStepperMPO::initUampo(const double U){
  ampo = Jampo;

  for(int i = 1; i <= sites.N(); ++i) {
    ampo += 0.5*U,"N(N-1)",i;
  }
}


void TimeStepperMPO::step(IQMPS& psi, const double from, const double to, const bool propagateForward){
  double U = 0.5*(from+to);
  initUampo(U);
  splitstep tOp;

  if (propagateForward) {
    auto expH1 = toExpH<IQTensor>(ampo,tstep*0.5*(1.0+Cplx_i));
    auto expH2 = toExpH<IQTensor>(ampo,tstep*0.5*(1.0-Cplx_i));
    tOp        = std::make_pair(expH2,expH1);
  }
  else {
    auto expH1 = toExpH<IQTensor>(ampo,-tstep*0.5*(1.0+Cplx_i));
    auto expH2 = toExpH<IQTensor>(ampo,-tstep*0.5*(1.0-Cplx_i));
    tOp        = std::make_pair(expH2,expH1);
  }

  doStep(psi, tOp);
}

void TimeStepperMPO::doStep(IQMPS& psi, const splitstep& tOp){
  psi = exactApplyMPO(tOp.first,psi,args);
  psi = exactApplyMPO(tOp.second,psi,args);
  normalize(psi);
}
