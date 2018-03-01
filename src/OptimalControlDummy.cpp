#include "OptimalControlDummy.hpp"


OptimalControlDummy::OptimalControlDummy(vec& weights, double tstep)
  : weights(weights), tstep(tstep) {

}

double OptimalControlDummy::getCost(const vec& control){
  double cost = 0;

  for (size_t n = 0; n < weights.size(); n++) {
    double an = weights.at(n);
    double cn = 0;

    for (auto& ui : control){
      cn += pow(ui,n);
    }
    cost += an*cn;
  }
  return cost*tstep;
}


vecpair OptimalControlDummy::getAnalyticGradient(const vec& control){
  vec grad;
  grad.reserve(control.size());

  for (auto& ui : control){
    double gi = 0;
    for (size_t n = 1; n < weights.size(); n++) {
      gi += weights.at(n)*n*pow(ui,n-1)*tstep;
    }
    grad.push_back(gi);
  }
  double cost = getCost(control);

  return std::make_pair(cost,grad);
}


vecpair OptimalControlDummy::getNumericGradient(const vec& control){
  auto newControl = control;
  double Jp, Jm;
  double epsilon = 1e-5;
  std::vector<double> g;
  g.reserve(control.size());

  for (auto& ui : newControl){
    ui        += epsilon;
    Jp         = getCost(newControl);

    ui        -= 2.0*epsilon;
    Jm         = getCost(newControl);

    ui        += epsilon;
    g.push_back((Jp-Jm)/(2.0*epsilon));
  }
  double cost = getCost(control);

  return std::make_pair(cost,g);
}


double OptimalControlDummy::getCost(const ControlBasis& bControl){
  return getCost(bControl.convControl() );
}


vecpair OptimalControlDummy::getAnalyticGradient(const ControlBasis& bControl){
  auto result = getAnalyticGradient(bControl.convControl());
  return std::make_pair(result.first,bControl.convGrad(result.second));
}


vecpair OptimalControlDummy::getNumericGradient(const ControlBasis& bControl){
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
