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

  for (auto &a : ampo) {
    auto expH1 = toExpH<IQTensor>(a,tau*0.5*(1.0+Cplx_i));
    auto expH2 = toExpH<IQTensor>(a,tau*0.5*(1.0-Cplx_i));

    psi = exactApplyMPO(expH2,psi,args);
    psi = exactApplyMPO(expH1,psi,args);
    normalize(psi);
    // fitApplyMPO(psi,expH2,psi,args);
    // fitApplyMPO(psi,expH1,psi,args);
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
  psi_t.reserve(nt+1);
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

template <class Iterable, class Tensor>
std::vector<IQMPS>
gateTEvolve2(Iterable const& gatelistlist,
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
  psi_t.reserve(nt+1);
  Real tsofar = 0;
  Real tot_norm = psi.normalize();
  psi_t.push_back(psi);
  if(verbose)
      {
      printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
      }
  psi.position(gatelistlist.front().front().i1());
  for(int tt = 1; tt <= nt; ++tt)
      {
      auto gatelist = gatelistlist[tt-1];
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
TimeEvolve2(const Iterable& gatelist,
          Real ttotal,
          Real tstep,
          MPSt<Tensor>& psi,
          const Args& args)
    {
    TEvolObserver obs(args);
    return gateTEvolve2(gatelist,ttotal,tstep,psi,obs,args);
}

template <class Iterable, class Tensor>
std::vector<IQMPS>
mixedTEvolve(Iterable const& gatelist,
          std::vector< MPOt<IQTensor> > const& MPOlist,
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
  psi_t.reserve(nt+1);
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
      bool forward = true;
      while(g != gatelist.end())
          {
          auto i1 = g->i1();
          auto i2 = g->i2();
          auto AA = psi.A(i1)*psi.A(i2)*g->gate();
          AA.mapprime(1,0,Site);

          if ( i1%2 && forward)
              {
              AA *= (MPOlist[tt-1]).A(i1)*(MPOlist[tt-1]).A(i2);
              AA.mapprime(1,0,Site);
              }

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
                  // std::cout << i1 << '\n';
                  psi.svdBond(i1,AA,Fromleft,args);
                  // auto AU = psi.A(i1)*(MPOlist[tt-1]).A(i1);
                  // AU.mapprime(1,0,Site);
                  // psi.setA(i1,AU);
                  psi.position(ni1); //does no work if position already ni1
                  }
              else
                  {
                  forward = false;
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
TimeEvolveMixed(const Iterable& gatelist,
          const std::vector< MPOt<IQTensor> >& MPOlist,
          Real ttotal,
          Real tstep,
          MPSt<Tensor>& psi,
          const Args& args)
    {
    TEvolObserver obs(args);
    return mixedTEvolve(gatelist,MPOlist,ttotal,tstep,psi,obs,args);
}

template <class Iterable, class Tensor>
std::vector<IQMPS>
mixedTEvolve2(Iterable const& gatelist,
          std::vector< std::vector<IQTensor> > const& expUlist,
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
  psi_t.reserve(nt+1);
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
      bool forward = true;
      while(g != gatelist.end())
          {
          auto i1 = g->i1();
          auto i2 = g->i2();
          auto AA = psi.A(i1)*psi.A(i2)*g->gate();
          AA.mapprime(1,0,Site);

          // if ( i1%2 && forward)
          //     {
          //     AA *= (MPOlist[tt-1]).A(i1)*(MPOlist[tt-1]).A(i2);
          //     AA.mapprime(1,0,Site);
          //     }

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
                  auto AU = psi.A(i1)*expUlist[tt-1][i1-1];
                  AU.mapprime(1,0,Site);
                  psi.setA(i1,AU);
                  psi.position(ni1); //does no work if position already ni1
                  }
              if(ni1 < i2 && !forward)
                  {
                  psi.svdBond(i1,AA,Fromright,args);
                  psi.position(ni2); //does no work if position already ni2
                  }
              if(ni1 < i2 && forward) {
                forward = false;
                psi.svdBond(i1,AA,Fromright,args);
                auto AU1 = psi.A(i1)*expUlist[tt-1][i1-1];
                auto AU2 = psi.A(i2)*expUlist[tt-1][i2-1];
                AU1.mapprime(1,0,Site);
                AU2.mapprime(1,0,Site);
                psi.setA(i1,AU1);
                psi.setA(i2,AU2);
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
TimeEvolveMixed2(const Iterable& gatelist,
          const std::vector< std::vector<IQTensor> >& expUlist,
          Real ttotal,
          Real tstep,
          MPSt<Tensor>& psi,
          const Args& args)
    {
    TEvolObserver obs(args);
    return mixedTEvolve2(gatelist,expUlist,ttotal,tstep,psi,obs,args);
}

template <class Iterable, class Tensor>
std::vector<IQMPS>
mixedTEvolve3(Iterable const& gatelist,
          std::vector< MPOt<IQTensor> > const& expUlist,
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
  psi_t.reserve(nt+1);
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
      bool forward = true;
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
                  auto AU = psi.A(i1)*(expUlist[tt-1]).A(i1);
                  AU.mapprime(1,0,Site);
                  psi.setA(i1,AU);
                  psi.position(ni1); //does no work if position already ni1
                  }
              if (ni1 < i2 && !forward)
                  {
                  psi.svdBond(i1,AA,Fromright,args);
                  psi.position(ni2); //does no work if position already ni2
                  }
              if (ni1 < i2 && forward)
                  {
                  forward = false;
                  psi.svdBond(i1,AA,Fromright,args);
                  auto AU1 = psi.A(i1)*(expUlist[tt-1]).A(i1);
                  auto AU2 = psi.A(i2)*(expUlist[tt-1]).A(i2);
                  AU1.mapprime(1,0,Site);
                  AU2.mapprime(1,0,Site);
                  psi.setA(i1,AU1);
                  psi.setA(i2,AU2);
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
TimeEvolveMixed3(const Iterable& gatelist,
          const std::vector< MPOt<IQTensor> >& expUlist,
          Real ttotal,
          Real tstep,
          MPSt<Tensor>& psi,
          const Args& args)
    {
    TEvolObserver obs(args);
    return mixedTEvolve3(gatelist,expUlist,ttotal,tstep,psi,obs,args);
}

} // end namespace
#endif
