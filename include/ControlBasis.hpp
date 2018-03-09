#ifndef CONTROLBASIS_H
#define CONTROLBASIS_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <armadillo>

typedef std::vector<double> stdvec;
// using mat = std::vector<std::vector<double> >;

class ControlBasis{
  // control has form u(t_i) = u0(t_i) + S(t_i)*sum(c_1*f_1(t_i) + ... + c_M*f_M(t_i))
private:

  arma::vec u0;
  arma::vec S;
  arma::vec c;
  arma::mat f;
  arma::mat ft;
  double dt;
  size_t M, N;


public:
  ControlBasis(arma::vec& u0, arma::vec& S, arma::mat& f, double dt);

  size_t getM() const;
  size_t getN() const;
  double getFij(size_t i, size_t j) const;
  void fmat2array(double* array);
  stdvec getCArray() const;
  void setCArray(const stdvec& cVec);
  void setCArray(const double* cArray, size_t size);
  stdvec convControl() const;
  void convControl(double* u, size_t size);
  stdvec convGrad(const stdvec& gradu) const;

};

#endif
