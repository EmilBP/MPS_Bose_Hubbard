#include "TimeStepperTEBDfast.hpp"

TimeStepperTEBDfast::TimeStepperTEBDfast(const SiteSet& sites, const double J, const double tstep, const Args& args)
  : J(J), sites(sites), args(args) {
  setTstep(tstep);
}

void TimeStepperTEBDfast::initJGates(const double J){
  using Gate = BondGate<IQTensor>;

  JGates_forwards.clear();
  JGates_backwards.clear();

  for(int i = 1; i < sites.N(); ++i)
      {
      auto hterm = -J*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -J*sites.op("Adag",i)*sites.op("A",i+1);

      auto gf = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      auto gb = Gate(sites,i,i+1,Gate::tReal,-tstep/2.,hterm);
      JGates_forwards.push_back(gf);
      JGates_backwards.push_back(gb);
      }
  for(int i = sites.N()-1; i >= 1; --i)
      {
      auto hterm = -J*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -J*sites.op("Adag",i)*sites.op("A",i+1);

      auto gf = Gate(sites,i,i+1,Gate::tReal,tstep/2.,hterm);
      auto gb = Gate(sites,i,i+1,Gate::tReal,-tstep/2.,hterm);
      JGates_forwards.push_back(gf);
      JGates_backwards.push_back(gb);
      }
}


void TimeStepperTEBDfast::setTstep(const double tstep_){
  tstep = tstep_;
  initJGates(J);
}

double TimeStepperTEBDfast::getTstep(){
  return tstep;
}

void TimeStepperTEBDfast::initUGates(const double U){
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

void TimeStepperTEBDfast::step(IQMPS& psi, const double from, const double to, const bool propagateForward){
  double U = 0.5*(from+to);

  if (propagateForward) {
    initUGates(U);
    doStep(psi, JGates_forwards, UGates);
  }
  else {
    initUGates(-U);
    doStep(psi, JGates_backwards, UGates);
  }
}

void TimeStepperTEBDfast::doStep(IQMPS& psi, const GateList& JGates, const std::vector<IQTensor> UGates){
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
              auto i1index = sites.si(i1);
              IQTensor A(i1index) ,B;
              denmatDecomp(AA,A,B,Fromleft,args);
              auto AU = A*UGates[i1-1];
              AU.mapprime(1,0,Site);
              B.mapprime(1,0,Site);
              psi.setA(i1,AU);
              psi.setA(i2,B);
              psi.position(ni1); //does no work if position already ni1
              }
          if(ni1 < i2 && !forward)
              {
              auto i1index = sites.si(i1);
              IQTensor A(i1index) ,B;
              denmatDecomp(AA,A,B,Fromright,args); //check if A and B need to be switched
              A.mapprime(1,0,Site);
              B.mapprime(1,0,Site);
              psi.setA(i1,A);
              psi.setA(i2,B);
              psi.position(ni2); //does no work if position already ni2
              }
          if(ni1 < i2 && forward) {
            forward = false;
            auto i1index = sites.si(i1);
            IQTensor A(i1index) ,B;
            denmatDecomp(AA,A,B,Fromright,args);
            auto AU1 = A*UGates[i1-1];
            auto AU2 = B*UGates[i2-1];
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
