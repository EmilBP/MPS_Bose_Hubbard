#ifndef CONTROLBASIS_H
#define CONTROLBASIS_H

#include <vector>
#include <iterator>
#include <algorithm>


using vec = std::vector<double>;
using mat = std::vector<std::vector<double> >;

class ControlBasis{
  // control has form u(t_i) = u0(t_i) + S(t_i)*sum(c_1*f_1(t_i) + ... + c_M*f_M(t_i))
private:

  vec u0;
  vec S;
  vec c;
  mat f;
  double dt;
  size_t M, N, controlIndex;


public:
  ControlBasis(vec& u0, vec& S, mat& f, double dt);

  size_t getM() const;
  size_t getN() const;
  size_t getControlIndex() const;
  double getFij(size_t i, size_t j) const;
  vec getCArray() const;
  void setCArray(const vec& cArray);
  vec convControl() const;
  vec convGrad(const vec& gradu) const;

  void exportParameters(vec& u_, vec& u0_, vec& S_, vec& c_, mat& f_);
};

#endif
