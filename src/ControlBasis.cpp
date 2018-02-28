#include "ControlBasis.hpp"

ControlBasis::ControlBasis(vec& u0, vec& S, mat& f)
  : u0(u0), S(S), f(f) {

}

vec ControlBasis::operator()(const vec& c){
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
