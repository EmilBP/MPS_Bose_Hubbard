#include "OptimalControl.hpp"
#include "TimeStepperTEBD.hpp"

template<class TimeStepper>
OptimalControl<TimeStepper>::OptimalControl(IQMPS& psi_target, IQMPS& psi_init, TimeStepper& timeStepper, double gamma, double tstep)
  : psi_target(psi_target), psi_init(psi_init), gamma(gamma), tstep(tstep), timeStepper(timeStepper) {

}

template<class TimeStepper>
double OptimalControl<TimeStepper>::getFidelity(const std::vector<double>& control){
  auto psi0 = psi_init;

  for (size_t i = 0; i < control.size()-1; i++) {
    timeStepper.step(psi0,control.at(i),control.at(i+1));
  }

  double re, im;
  overlap(psi_target,psi0,re,im);
  return (re*re+im*im);
}

template<class TimeStepper>
double OptimalControl<TimeStepper>::getRegularisation(const std::vector<double>& control){

  double tmp = 0;
  for (size_t i = 0; i < control.size()-1; i++) {

    double diff = control.at(i+1)-control.at(i);
    tmp += diff*diff/tstep;
  }

  return gamma/2.0*tmp;
}

template<class TimeStepper>
double OptimalControl<TimeStepper>::getCost(const std::vector<double>& control){
  return getRegularisation(control);

  // return 0.5*(1-getFidelity(control)) + getRegularisation(control);
}

template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getRegGrad(const std::vector<double>& control){
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
std::vector<double> OptimalControl<TimeStepper>::getFidelityGrad(const std::vector<double>& control){
  std::vector<double> g;
  // g.reserve(control.size());
  //
  // auto i1 = begin(chi_t), i2 = begin(psi_t);
  // auto i3 = begin(control), e = end(control);
  //
  // while (i3 != e) {
  //   g.emplace_back(tstep*overlapC(*i1,IQMPO(MPOstrat.dHdu(*i3)),*i2).real() );
  //   ++i1;
  //   ++i2;
  //   ++i3;
  // }
  return g;
}

template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getAnalyticGradient(const std::vector<double>& control){
  return getRegGrad(control);
}


template<class TimeStepper>
std::vector<double> OptimalControl<TimeStepper>::getNumericGradient(const std::vector<double>& control){

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
//
// std::vector<double> OptimalControl::getAnalyticGradientTest(std::vector<double>& control, double dt, const Args& args){
//   auto ampo = updateHamiltonian(control);
//   psi_t = TimeEvolve(psi_init, ampo, dt*Cplx_i, args);
//   auto chi_T = -Cplx_i * psi_target*overlapC(psi_target,psi_t.back());
//   chi_t = TimeEvolveBack(chi_T, ampo, dt*Cplx_i, args);
//
//   return getAnalyticGradient(control,dt);
// }
//
// std::vector<double> OptimalControl::getNumericGradient(
//       std::vector<double>& control, double epsilon,
//       double dt, const Args& args)
// {
//   double Jp, Jm;
//   std::vector<double> g;
//   std::vector<IQMPS> psi_temp;
//   std::vector<AutoMPO> tempMPO;
//   g.reserve(control.size());
//
//   for (auto& ui : control){
//     ui        += epsilon;
//     tempMPO    = updateHamiltonian(control);
//     psi_temp   = TimeEvolve(psi_init,tempMPO,dt*Cplx_i,args);
//     Jp         = getCost(psi_temp.back(),control);
//
//     ui        -= 2*epsilon;
//     tempMPO    = updateHamiltonian(control);
//     psi_temp   = TimeEvolve(psi_init,tempMPO,dt*Cplx_i,args);
//     Jm         = getCost(psi_temp.back(),control);
//
//     ui += epsilon;
//
//     g.emplace_back((Jp-Jm)/(2.0*epsilon));
//   }
//
//   return g;
// }
