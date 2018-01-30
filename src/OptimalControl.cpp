#include "OptimalControl.hpp"
#include "TimeStepperTEBD.hpp"

template<class TimeStepper>
OptimalControl<TimeStepper>::OptimalControl(IQMPS& psi_target, IQMPS& psi_init, TimeStepper& timeStepper, MPOt<IQTensor>& dHdU, double gamma, double tstep)
  : psi_target(psi_target), psi_init(psi_init), gamma(gamma), tstep(tstep), timeStepper(timeStepper), dHdU(dHdU) {

}

template<class TimeStepper>
double OptimalControl<TimeStepper>::getFidelity(vec& control){
  auto psi0 = psi_init;

  for (size_t i = 0; i < control.size()-1; i++) {
    timeStepper.step(psi0,control.at(i),control.at(i+1));
  }

  double re, im;
  overlap(psi_target,psi0,re,im);
  return (re*re+im*im);
}

template<class TimeStepper>
double OptimalControl<TimeStepper>::getRegularisation(vec& control){

  double tmp = 0;
  for (size_t i = 0; i < control.size()-1; i++) {

    double diff = control.at(i+1)-control.at(i);
    tmp += diff*diff/tstep;
  }

  return gamma/2.0*tmp;
}

template<class TimeStepper>
double OptimalControl<TimeStepper>::getCost(vec& control){

  return 0.5*(1-getFidelity(control)) + getRegularisation(control);
}

template<class TimeStepper>
void OptimalControl<TimeStepper>::calcPsi(vec& control){
  auto psi0 = psi_init;
  psi_t.clear();
  psi_t.push_back(psi0);

  for (size_t i = 0; i < control.size()-1; i++) {
    timeStepper.step(psi0,control.at(i),control.at(i+1));
    psi_t.push_back(psi0);
  }
}

template<class TimeStepper>
void OptimalControl<TimeStepper>::calcChi(vec& control){
  auto chiT = psi_target;
  // auto chiT = -Cplx_i * psi_target*overlapC(psi_target,psi_t.back());
  chi_t.clear();
  chi_t.push_back(chiT);

  for (size_t i = control.size()-1; i > 0; i--) {
    timeStepper.step(chiT,control.at(i),control.at(i-1));
    chi_t.push_back(chiT);
  }

  std::reverse(chi_t.begin(),chi_t.end());
}

template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getRegGrad(vec& control){
  std::vector<double> del;
  del.reserve(control.size());

  del.push_back(-gamma*(-5.0*control.at(1) + 4.0*control.at(2) - control.at(3)
                  + 2.0*control.at(0))/tstep/tstep);

  for (size_t i = 1; i < control.size()-1; i++) {
    del.push_back(-gamma*(control.at(i+1) + control.at(i-1) - 2.0*control.at(i))/tstep/tstep);
  }

  del.push_back( -gamma*(-5.0*control.at(control.size()-2) + 4.0*control.at(control.size()-3)
                  - control.at(control.size()-4) + 2.0*control.at(control.size()-1))/tstep/tstep);

  return del;
}

template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getFidelityGrad(vec& control){
  std::vector<double> g;
  g.reserve(control.size());

  calcPsi(control);
  calcChi(control);

  for (size_t i = 0; i < control.size(); i++) {
    g.push_back(-tstep*overlapC(chi_t.at(i),dHdU,psi_t.at(i)).real() );
  }

  return g;
}

template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getAnalyticGradient(vec& control){

  auto FGrad = getFidelityGrad(control);
  auto RGrad = getRegGrad(control);

  for (size_t i = 0; i < FGrad.size(); i++) {
    FGrad.at(i) += RGrad.at(i);
  }
  return FGrad;
  // return getFidelityGrad(control) + getRegGrad(control);
}


template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getNumericGradient(vec& control){

  auto newControl = control;
  double Jp, Jm;
  double epsilon = 1e-6;
  std::vector<double> g;
  g.reserve(control.size());

  for (auto& ui : newControl){
    ui        += epsilon;
    Jp         = getCost(newControl);

    ui        -= 2*epsilon;
    Jm         = getCost(newControl);

    ui        += epsilon;

    g.push_back((Jp-Jm)/(2.0*epsilon));
  }

  return g;
}

template class OptimalControl<TimeStepperTEBD>;
