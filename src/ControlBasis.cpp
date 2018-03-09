#include "ControlBasis.hpp"

ControlBasis::ControlBasis(arma::vec& u0, arma::vec& S, arma::mat& f, double dt)
  : u0(u0), S(S), f(f), ft(f.t()), dt(dt), M(f.n_cols), N(u0.size()) {

  c = arma::zeros<arma::vec>(M);
}

stdvec ControlBasis::getCArray() const{
  // return as std::vec
  return arma::conv_to< stdvec >::from(c);
}

size_t ControlBasis::getM() const{
  return M;
}

size_t ControlBasis::getN() const{
  return N;
}

double ControlBasis::getFij(size_t i, size_t j) const{
  return f(i,j);
}

void ControlBasis::fmat2array(double* array) {
  // return f as a single, long array double*
  // uses transpose as arma::mat stored as columns, but Ipopt takes rows
  array = ft.memptr();
}


void ControlBasis::setCArray(const stdvec& cVec){
  // convert std::vec input to arma::vec member
  c = arma::conv_to< arma::vec >::from(cVec);
}

void ControlBasis::setCArray(const double* cArray, size_t size) {
  // convert double* input to arma::vec member
  c = arma::vec(cArray, size);
}

stdvec ControlBasis::convControl() const{
  // u = u0 + S*sum(fn*cn)
  // convert output to std::vec
  return arma::conv_to< stdvec >::from( u0+S%(f*c) );
}

void ControlBasis::convControl(double* u, size_t size) {
  // calculate  arma::vec u and return as double*
  arma::vec uv = u0+S%(f*c);
  u = uv.memptr();
}


stdvec ControlBasis::convGrad(const stdvec& gradu) const{
  // convert std::vec input to arma::vec for faster calculations
  // convert back to std::vec and return
  arma::vec G = arma::conv_to< arma::vec >::from(gradu);
  return arma::conv_to< stdvec >::from( ft*(G%S) );
}
