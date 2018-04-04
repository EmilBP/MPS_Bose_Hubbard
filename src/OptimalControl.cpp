#include "OptimalControl.hpp"
#include "TimeStepperTEBD.hpp"
#include "TimeStepperTEBDfast.hpp"
#include "TimeStepperTEBDnew.hpp"
#include "HamiltonianBH.hpp"

template<class TimeStepper, class Hamiltonian>
OptimalControl<TimeStepper,Hamiltonian>::OptimalControl(IQMPS& psi_target, IQMPS& psi_init, TimeStepper& timeStepper, Hamiltonian& hamil, double gamma)
  : psi_target(psi_target), psi_init(psi_init), gamma(gamma), timeStepper(timeStepper), hamil(hamil), tstep(timeStepper.getTstep()) {

}

template<class TimeStepper, class Hamiltonian>
double OptimalControl<TimeStepper,Hamiltonian>::getFidelity(const vec& control){
  auto psi0 = psi_init;

  for (size_t i = 0; i < control.size()-1; i++) {
    timeStepper.step(psi0,control.at(i),control.at(i+1));
  }

  double re, im;
  overlap(psi_target,psi0,re,im);
  return (re*re+im*im);
}

template<class TimeStepper, class Hamiltonian>
double OptimalControl<TimeStepper,Hamiltonian>::getRegularisation(const vec& control){

  double tmp = 0;
  for (size_t i = 0; i < control.size()-1; i++) {

    double diff = control.at(i+1)-control.at(i);
    tmp += diff*diff/tstep;
  }

  return gamma/2.0*tmp;
}

template<class TimeStepper, class Hamiltonian>
double OptimalControl<TimeStepper,Hamiltonian>::getCost(const vec& control){

  return 0.5*(1-getFidelity(control)) + getRegularisation(control);
}

template<class TimeStepper, class Hamiltonian>
void OptimalControl<TimeStepper,Hamiltonian>::calcPsi(const vec& control){

  const bool propagateForward = true;
  auto psi0 = psi_init;
  psi_t.clear();
  psi_t.push_back(psi0);

  for (size_t i = 0; i < control.size()-1; i++) {
    timeStepper.step(psi0,control.at(i),control.at(i+1),propagateForward);
    psi_t.push_back(psi0);
  }
}

template<class TimeStepper, class Hamiltonian>
void OptimalControl<TimeStepper,Hamiltonian>::calcChi(const vec& control){

  const bool propagateForward = false;
  auto chiT = Cplx_i*psi_target;
  chi_t.clear();
  chi_t.push_back(chiT);

  for (size_t i = control.size()-1; i > 0; i--) {
    timeStepper.step(chiT,control.at(i),control.at(i-1),propagateForward);
    chi_t.push_back(chiT);
  }

  std::reverse(chi_t.begin(),chi_t.end());
}

template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::getRegPlusRegGrad(const vec& control){

  auto reg = getRegularisation(control);
  std::vector<double> del;
  del.reserve(control.size());

  del.push_back(-gamma*(-5.0*control.at(1) + 4.0*control.at(2) - control.at(3)
                  + 2.0*control.at(0))/tstep/tstep);

  for (size_t i = 1; i < control.size()-1; i++) {
    del.push_back(-gamma*(control.at(i+1) + control.at(i-1) - 2.0*control.at(i))/tstep/tstep);
  }

  del.push_back( -gamma*(-5.0*control.at(control.size()-2) + 4.0*control.at(control.size()-3)
                  - control.at(control.size()-4) + 2.0*control.at(control.size()-1))/tstep/tstep);

  return std::make_pair(reg,del);
}

template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::getFidelityPlusFidelityGrad(const vec& control){


  std::vector<double> g;
  g.reserve(control.size());

  calcPsi(control);
  calcChi(control);

  auto overlapFactor = overlapC(psi_t.back(),psi_target);

  for (size_t i = 0; i < control.size(); i++) {
    g.push_back( -(overlapC( chi_t.at(i) , hamil.dHdU(control.at(i)) , psi_t.at(i) )*overlapFactor ).real() );
  }

  double re, im;
  overlap(psi_target,psi_t.back(),re,im);
  double fidelity = (re*re+im*im);

  return std::make_pair(fidelity,g);
}

template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::getAnalyticGradient(const vec& control){

  auto FGrad = getFidelityPlusFidelityGrad(control);
  auto RGrad = getRegPlusRegGrad(control);

  for (size_t i = 0; i < FGrad.second.size(); i++) {
    FGrad.second.at(i) += RGrad.second.at(i);
  }

  double cost = 0.5*(1-FGrad.first) + RGrad.first;
  FGrad.first = cost;

  return FGrad;
}


template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::getNumericGradient(const vec& control){

  auto newControl = control;
  double Jp, Jm;
  double epsilon = 1e-5;
  std::vector<double> g;
  g.reserve(control.size());

  size_t count = 0;

  for (auto& ui : newControl){
    ui        += epsilon;
    Jp         = getCost(newControl);

    ui        -= 2.0*epsilon;
    Jm         = getCost(newControl);

    ui        += epsilon;

    std::cout << count++ << '\n';

    g.push_back((Jp-Jm)/(2.0*epsilon));
  }
  double cost = getCost(control);

  return std::make_pair(cost,g);
}

template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::checkCostPlusFidelity(const vec& control){

  std::vector<double> fid;
  fid.reserve(control.size());

  calcPsi(control);
  calcChi(control);

  double re, im;
  for (size_t i = 0; i < control.size(); i++) {
    overlap(chi_t.at(i),psi_t.at(i),re,im);
    fid.push_back(re*re+im*im);
  }

  overlap(psi_target,psi_t.back(),re,im);
  double cost = 0.5*(1-(re*re+im*im)) + getRegularisation(control);

  return std::make_pair(cost,fid);
}

template<class TimeStepper, class Hamiltonian>
double OptimalControl<TimeStepper,Hamiltonian>::getCost(const ControlBasis& bControl){

  return getCost(bControl.convControl());
}

template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::getAnalyticGradient(const ControlBasis& bControl){

  auto result = getAnalyticGradient(bControl.convControl());
  return std::make_pair(result.first,bControl.convGrad(result.second));
}

template<class TimeStepper, class Hamiltonian>
vecpair OptimalControl<TimeStepper,Hamiltonian>::getNumericGradient(const ControlBasis& bControl){
  auto newbControl = bControl;
  auto cArray = newbControl.getCArray();
  double Jp, Jm;
  double epsilon = 1e-5;
  std::vector<double> g;

  for (auto& ci : cArray){
    ci        += epsilon;
    newbControl.setCArray(cArray);
    Jp         = getCost(newbControl);

    ci        -= 2.0*epsilon;
    newbControl.setCArray(cArray);
    Jm         = getCost(newbControl);

    ci        += epsilon;
    newbControl.setCArray(cArray);
    g.push_back((Jp-Jm)/(2.0*epsilon));
  }
  double cost = getCost(newbControl);

  return std::make_pair(cost,g);
}

template<class TimeStepper, class Hamiltonian>
vec OptimalControl<TimeStepper,Hamiltonian>::getFidelityForAllT(const vec& control){

  std::vector<double> fid;
  fid.reserve(control.size());
  calcPsi(control);

  double re, im;
  for (size_t i = 0; i < control.size(); i++) {
    overlap(psi_target,psi_t.at(i),re,im);
    fid.push_back(re*re+im*im);
  }

  return fid;
}

template<class TimeStepper, class Hamiltonian>
vec OptimalControl<TimeStepper,Hamiltonian>::getFidelityForAllT(const ControlBasis& bControl){
  return getFidelityForAllT(bControl.convControl());
}

template<class TimeStepper, class Hamiltonian>
std::vector<IQMPS> OptimalControl<TimeStepper,Hamiltonian>::getPsit() const{
  return psi_t;
}

template class OptimalControl<TimeStepperTEBD,HamiltonianBH>;
template class OptimalControl<TimeStepperTEBDfast,HamiltonianBH>;
template class OptimalControl<TimeStepperTEBDnew,HamiltonianBH>;
