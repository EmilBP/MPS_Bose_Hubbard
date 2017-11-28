#ifndef TIMEEVOLVE_H
#define TIMEEVOLVE_H

#include "itensor/all.h"
#include <vector>
#include <iterator>
#include <algorithm>


namespace itensor{

//
//  time evolution methods
//

inline std::vector<IQMPS> TimeEvolve(IQMPS& psi, std::vector<AutoMPO>& ampo, Complex tau, const Args& args){
  std::vector<IQMPS> psi_t;
  psi_t.reserve(ampo.size()+1);
  psi_t.push_back(psi);

  for (auto &a : ampo) {
    auto expH1 = toExpH<IQTensor>(a,tau*0.5*(1.0+Cplx_i));
    auto expH2 = toExpH<IQTensor>(a,tau*0.5*(1.0-Cplx_i));

    fitApplyMPO(psi,expH2,psi,args);
    fitApplyMPO(psi,expH1,psi,args);
    psi_t.push_back(psi);
  }

  // psi_t.pop_back();

  return psi_t;
}

inline std::vector<IQMPS> TimeEvolve(IQMPS& psi, AutoMPO& ampo, Complex tau, size_t Nsteps, const Args& args){
  auto expH1 = toExpH<IQTensor>(ampo,tau*0.5*(1.0+Cplx_i));
  auto expH2 = toExpH<IQTensor>(ampo,tau*0.5*(1.0-Cplx_i));

  std::vector<IQMPS> psi_t;
  psi_t.reserve(Nsteps);
  psi_t.push_back(psi);

  for (size_t t = 0; t < Nsteps; ++t) {
    fitApplyMPO(psi,expH2,psi,args);
    fitApplyMPO(psi,expH1,psi,args);
    psi_t.push_back(psi);
  }

  return psi_t;
}

inline std::vector<IQMPS> TimeEvolveBack(IQMPS& psi, std::vector<AutoMPO>& ampo, Complex tau, const Args& args){
  std::vector<IQMPS> psi_t;
  psi_t.reserve(ampo.size()+1);
  psi_t.push_back(psi);

  for (auto it = ampo.rbegin(); it != ampo.rend(); ++it){
    auto expH1 = toExpH<IQTensor>((*it),-tau*0.5*(1.0+Cplx_i));
    auto expH2 = toExpH<IQTensor>((*it),-tau*0.5*(1.0-Cplx_i));

    fitApplyMPO(psi,expH2,psi,args);
    fitApplyMPO(psi,expH1,psi,args);
    psi_t.push_back(psi);
  }

  // psi_t.pop_back();

  std::reverse(psi_t.begin(),psi_t.end());
  return psi_t;
}




} // end namespace
#endif
