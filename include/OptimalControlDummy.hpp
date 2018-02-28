#ifndef OPTIMALCONTROLDUMMY_H
#define OPTIMALCONTROLDUMMY_H

#include "IpIpoptApplication.hpp"
#include "itensor/all.h"
#include <vector>
#include <math.h>
#include <iterator>
#include <algorithm>



using vec = std::vector<double>;
using vecpair = std::pair<double, vec>;

class OptimalControlDummy{
private:
  vec weights;
  double tstep;



public:
  OptimalControlDummy(vec& weights, double tstep);

  double getCost(const vec& control);
  vecpair getAnalyticGradient(const vec& control);
  vecpair getNumericGradient(const vec& control);

};

#endif
