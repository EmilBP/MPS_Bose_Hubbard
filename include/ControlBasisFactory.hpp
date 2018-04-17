#ifndef CONTROLBASISFACTORY_HPP
#define CONTROLBASISFACTORY_HPP

#include "ControlBasis.hpp"
#include <vector>
#include <assert.h>
#include <armadillo>
#include "math.h"

#define PI 3.14159265
typedef std::vector<double> stdvec;

class ControlBasisFactory{

private:
  static arma::vec sigmoid(arma::vec& x, double k, double offset);

public:

  static ControlBasis buildCBsin(stdvec& u0, double tstep, double T, size_t M);

};

arma::vec ControlBasisFactory::sigmoid(arma::vec& x, double k, double offset){
  arma::vec S = x;

  for (size_t i = 0; i < S.n_rows; i++) {
    S(i) = 1.0/(1+exp(-k*(x(i)-offset) ) );
  }

  return S;
}

ControlBasis ControlBasisFactory::buildCBsin(stdvec& u0, double tstep, double T, size_t M){
  size_t N = u0.size();
  assert( N-(1 + T/tstep) < 1e-5 );
  arma::vec x    = arma::linspace<arma::vec>(0,100,N);
  arma::vec S    = sigmoid(x,8.0,1.1);
  arma::vec S2   = sigmoid(x,-8.0,100-1.1);

  for (size_t i = S.n_rows/2; i < S.n_rows; i++) {
    S(i) = S2(i);
  }
  S(0)          = 0;
  S(S.n_rows-1) = 0;

  arma::mat f(N,M);
  for (size_t i = 0; i < N; i++) {
    for (size_t n = 1; n <= M; n++) {
      f(i,n-1) = sin(n*PI*tstep*i/T);
    }
  }

  arma::vec u0a = arma::conv_to< arma::vec >::from(u0);
  return ControlBasis(u0a,S,f,tstep);
}


#endif
