#include "ControlBasis.hpp"

ControlBasis::ControlBasis(arma::vec& u0, arma::vec& S, arma::mat& f, double dt)
  : u0(u0), S(S), f(f), ft(f.t()), dt(dt), M(f.n_cols), N(u0.size()) {

  constraints = (arma::diagmat(S)*f).t();
  c           = arma::zeros<arma::vec>(M);
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

stdvec ControlBasis::getU0() const{
  // return as std::vec
  return arma::conv_to< stdvec >::from(u0);
}

double ControlBasis::getFij(size_t i, size_t j) const{
  return f(i,j);
}

void ControlBasis::getConstraintJacobian(double* array) {
  // return f as a single, long array double*
  // uses transpose as arma::mat stored as columns, but Ipopt takes rows
  double* farray = constraints.memptr();
  std::copy(farray, farray+N*M, array);
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

void ControlBasis::convControl(double* u) {
  // calculate  arma::vec u and return as double*
  arma::vec uv = u0+S%(f*c);
  double* up = uv.memptr();
  std::copy(up, up+N, u);
}


stdvec ControlBasis::convGrad(const stdvec& gradu) const{
  // convert std::vec input to arma::vec for faster calculations
  // convert back to std::vec and return
  arma::vec G = arma::conv_to< arma::vec >::from(gradu);
  return arma::conv_to< stdvec >::from( ft*(G%S) );
}


void ControlBasis::exportParameters(){
  arma::mat vecdata = arma::join_horiz( u0+S%(f*c) , u0);
  vecdata = arma::join_horiz(vecdata , S);

  arma::mat matdata = arma::join_vert(f, c.t());

  std::string filename1 = "CBsinVecdata_M" + std::to_string(M) + ".txt";
  std::string filename2 = "CBsinMatdata_M" + std::to_string(M) + ".txt";

  vecdata.save(filename1,arma::raw_ascii);
  matdata.save(filename2,arma::raw_ascii);
}
