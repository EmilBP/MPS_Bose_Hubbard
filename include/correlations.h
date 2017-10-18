#ifndef CORRELATIONS_H
#define CORRELATIONS_H

#include "itensor/all.h"

using namespace itensor;

class correlations{
public:
  static Cplx correlationFunction(SiteSet& sites, IQMPS& psi, std::string const& opname1, int i, std::string const& opname2, int j);
  static Real correlationTerm(SiteSet sites, IQMPS psi, std::string const& opname1, std::string const& opname2);
};

#endif
