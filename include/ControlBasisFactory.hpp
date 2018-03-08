#ifndef CONTROLBASISFACTORY_HPP
#define CONTROLBASISFACTORY_HPP

#include "ControlBasis.hpp"
#include <vector>
#include <assert.h>
#include "math.h"

#define PI 3.14159265

class ControlBasisFactory{

private:
  static std::vector<double> linspace(double a, double b, int n);
  static std::vector<double> sigmoid(std::vector<double>& x, double k, double offset);

public:

  static ControlBasis buildCBsin(std::vector<double>& u0, double tstep, double T, size_t M);

};

std::vector<double> ControlBasisFactory::linspace(double a, double b, int n){
  std::vector<double> array;
  double step = (b-a) / (n-1);

  while(a <= b + 1e-7) {
      array.push_back(a);
      a += step;           // could recode to better handle rounding errors
  }
  return array;
}

std::vector<double> ControlBasisFactory::sigmoid(std::vector<double>& x, double k, double offset){
  std::vector<double> S;
  for (auto& xval : x){
    S.push_back(1.0/(1+exp(-k*(xval-offset))));
  }
  return S;
}

ControlBasis ControlBasisFactory::buildCBsin(std::vector<double>& u0, double tstep, double T, size_t M){
  assert( u0.size()-(1 + T/tstep) < 1e-5 );
  auto x    = linspace(0,100,u0.size());
  auto S    = sigmoid(x,3.0,2.5);
  auto S2   = sigmoid(x,-3.0,100-2.5);

  for (size_t i = S.size()/2; i < S.size(); i++) {
    S.at(i) = S2.at(i);
  }
  S.front() = 0;
  S.back()  = 0;

  std::vector<std::vector<double> > f;
  for (size_t i = 0; i < u0.size(); i++) {
    std::vector<double> fi;
    for (size_t n = 1; n <= M; n++) {
      fi.push_back( sin(n*PI*tstep*i/T) );
    }
    f.push_back(fi);
  }

  return ControlBasis(u0,S,f,tstep);
}


#endif
