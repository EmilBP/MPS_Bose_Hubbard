#include "TimeStepperTEBD.hpp"

TimeStepperTEBD::TimeStepperTEBD(const SiteSet& sites, const double J, const double tstep, const Args& args)
  : J(J), sites(sites), args(args) {
  setTstep(tstep);
}

void TimeStepperTEBD::initJGates(const double J){
  using Gate = BondGate<IQTensor>;

  JGates.clear();
  for(int i = 1; i < sites.N(); ++i)
      {
      auto hterm = -J*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -J*sites.op("Adag",i)*sites.op("A",i+1);

      auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      JGates.push_back(g);
      }
  for(int i = sites.N()-1; i >= 1; --i)
      {
      auto hterm = -J*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -J*sites.op("Adag",i)*sites.op("A",i+1);

      auto g = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      JGates.push_back(g);
      }
}

void TimeStepperTEBD::setTstep(const double tstep_){
  tstep = tstep_;
  initJGates(J);
}

void TimeStepperTEBD::initUGates(const double U){
  UGates.clear();

  for (int k = 1; k <= sites.N(); ++k) {
    auto s    = sites.si(k);
    auto sP   = prime(s);
    int HD    = s.nblock();

    IQTensor T(dag(s),sP);

    for (size_t i = 0; i < HD; i++) {
      T.set(s(i+1),sP(i+1), std::exp( -0.5*U*tstep*Cplx_i*i*(i-1) ) );
    }

    UGates.push_back(T);
  }
}

void TimeStepperTEBD::step(IQMPS& psi, const double from, const double to){
  double U = 0.5*(from+to);
  initUGates(U);

  const bool normalize = args.getBool("Normalize",true);
  auto g = JGates.begin();
  bool forward = true;
  while(g != JGates.end())
      {
      auto i1 = g->i1();
      auto i2 = g->i2();
      auto AA = psi.A(i1)*psi.A(i2)*g->gate();
      AA.mapprime(1,0,Site);

      ++g;
      if(g != JGates.end())
          {
          //Look ahead to next gate position
          auto ni1 = g->i1();
          auto ni2 = g->i2();
          //SVD AA to restore MPS form
          //before applying current gate
          if(ni1 >= i2)
              {
              psi.svdBond(i1,AA,Fromleft,args);
              auto AU = psi.A(i1)*UGates[i1-1];
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
            auto AU1 = psi.A(i1)*UGates[i1-1];
            auto AU2 = psi.A(i2)*UGates[i2-1];
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
      psi.normalize();
      }
}
