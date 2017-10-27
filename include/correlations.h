#ifndef CORRELATIONS_H
#define CORRELATIONS_H

#include "itensor/all.h"
#include <iostream>

using namespace itensor;

class correlations{
public:
  static Cplx correlationFunction(SiteSet const& sites, IQMPS& psi, std::string const& opname1, int i, std::string const& opname2, int j);
  static ITensor correlationMatrix(SiteSet const& sites, IQMPS& psi, std::string const& opname1, std::string const& opname2);
  static Real correlationTerm(SiteSet const& sites, IQMPS& psi, std::string const& opname1, std::string const& opname2);
};

#endif
