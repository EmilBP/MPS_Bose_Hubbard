#ifndef CONTROLBASISFACTORY_HPP
#define CONTROLBASISFACTORY_HPP

#include "ControlBasis.hpp"
#include <vector>
#include "math.h"

#define PI 3.14159265

class ControlBasisFactory{

private:
  static std::vector<double> linspace(double a, double b, int n);

public:

  static ControlBasis buildCRAB(double cstart, double cend, double tstep, double T, size_t M);

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

ControlBasis ControlBasisFactory::buildCRAB(double cstart, double cend, double tstep, double T, size_t M){
  auto u0 = linspace(cstart,cend,T/tstep+1);
  std::vector<double> S(u0.size(),1.0);
  S.front() = 0;
  S.back() = 0;

  std::vector<std::vector<double> > f;
  for (size_t i = 0; i < u0.size(); i++) {
    std::vector<double> fi;
    for (size_t n = 0; n < M; n++) {
      fi.push_back( sin(n*PI*tstep*i/T) );
    }
    f.push_back(fi);
  }

  return ControlBasis(u0,S,f,tstep);
}


#endif
