#ifndef TIMEEVOLVE_H
#define TIMEEVOLVE_H

#include "itensor/all.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/bondgate.h"
#include "itensor/mps/TEvolObserver.h"
#include <vector>
#include <time.h>
#include <iterator>
#include <algorithm>


namespace itensor{

//
//  time evolution methods
//

inline std::vector<IQMPS> TimeEvolve(IQMPS psi, std::vector<AutoMPO>& ampo, Complex tau, const Args& args){
  std::vector<IQMPS> psi_t;
  psi_t.reserve(ampo.size()+1);
  psi_t.push_back(psi);
  size_t index = 0;

  // std::cout << "Time Evolution" << '\n';
  // std::cout << index << "/" << ampo.size()+1 << '\r';
  for (auto &a : ampo) {
    clock_t begin = clock();

    auto expH1 = toExpH<IQTensor>(a,tau*0.5*(1.0+Cplx_i));
    auto expH2 = toExpH<IQTensor>(a,tau*0.5*(1.0-Cplx_i));
    clock_t end = clock();
    std::cout << "Runtime for exponent = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";
    begin = clock();
    fitApplyMPO(psi,expH2,psi,args);
    fitApplyMPO(psi,expH1,psi,args);
    end = clock();
    std::cout << "Runtime for applying = " <<  double(end - begin) / CLOCKS_PER_SEC << "\n";
    psi_t.push_back(psi);
    // std::cout << index++ << "/" << ampo.size()+1 << '\r';
  }
  // std::cout << '\n';

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

template <class Iterable, class Tensor>
std::vector<IQMPS>
gateTEvolve(Iterable const& gatelist,
          Real ttotal,
          Real tstep,
          MPSt<Tensor>& psi,
          Observer& obs,
          Args args)
{
  const bool verbose = args.getBool("Verbose",false);
  const bool normalize = args.getBool("Normalize",true);

  const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
  if(fabs(nt*tstep-ttotal) > 1E-9)
      {
      Error("Timestep not commensurate with total time");
      }

  std::vector<IQMPS> psi_t;
  psi_t.reserve(ttotal/tstep+1);
  Real tsofar = 0;
  Real tot_norm = psi.normalize();
  psi_t.push_back(psi);
  if(verbose)
      {
      printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
      }
  psi.position(gatelist.front().i1());
  for(int tt = 1; tt <= nt; ++tt)
      {
      auto g = gatelist.begin();
      while(g != gatelist.end())
          {
          auto i1 = g->i1();
          auto i2 = g->i2();
          auto AA = psi.A(i1)*psi.A(i2)*g->gate();
          AA.mapprime(1,0,Site);

          ++g;
          if(g != gatelist.end())
              {
              //Look ahead to next gate position
              auto ni1 = g->i1();
              auto ni2 = g->i2();
              //SVD AA to restore MPS form
              //before applying current gate
              if(ni1 >= i2)
                  {
                  psi.svdBond(i1,AA,Fromleft,args);
                  psi.position(ni1); //does no work if position already ni1
                  }
              else
                  {
                  psi.svdBond(i1,AA,Fromright,args);
                  psi.position(ni2); //does no work if position already ni2
                  }
              }
          else
              {
              //No next gate to analyze, just restore MPS form
              psi.svdBond(i1,AA,Fromright,args);
              }
          }

      if(normalize)
          {
          tot_norm *= psi.normalize();
          }

      psi_t.push_back(psi);
      tsofar += tstep;

      args.add("TimeStepNum",tt);
      args.add("Time",tsofar);
      args.add("TotalTime",ttotal);
      obs.measure(args);
      }
  if(verbose)
      {
      printfln("\nTotal time evolved = %.5f\n",tsofar);
      }

  return psi_t;

}

template <class Iterable, class Tensor>
inline std::vector<IQMPS>
TimeEvolve(const Iterable& gatelist, 
          Real ttotal,
          Real tstep,
          MPSt<Tensor>& psi,
          const Args& args)
    {
    TEvolObserver obs(args);
    return gateTEvolve(gatelist,ttotal,tstep,psi,obs,args);
}


} // end namespace
#endif
