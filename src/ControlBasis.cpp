#include "ControlBasis.hpp"

ControlBasis::ControlBasis(vec& u0, vec& S, mat& f, double dt)
  : u0(u0), S(S), f(f), dt(dt), M(f.front().size()) {

  c = std::vector<double>(M,0);
}

vec ControlBasis::getCArray() const{
  return c;
}

size_t ControlBasis::getM() const{
  return M;
}

void ControlBasis::setCArray(const vec& cArray){
  c = cArray;
}

vec ControlBasis::convControl() const{
  vec u = u0;
  auto Si = S.begin();
  auto fi = f.begin();

  for (auto& ui : u){
    auto cn = c.begin();
    double sum = 0;

    for (auto& fn : (*fi) ){
      sum += fn*(*cn++);
    }

    ui += (*Si)*sum;

    fi++; Si++;
  }
  return u;
}

vec ControlBasis::convGrad(const vec& gradu) const{
  vec gc;
  gc.reserve(M);

  for (size_t n = 0; n < M; n++) {
    double gc_n = 0;

    for (size_t i = 0; i < gradu.size(); i++) {
      gc_n += gradu.at(i)*S.at(i)*(f.at(i)).at(n);
    }
    gc.push_back(gc_n);
  }
  return gc;

}
