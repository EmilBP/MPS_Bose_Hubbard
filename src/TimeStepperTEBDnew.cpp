#include "TimeStepperTEBDnew.hpp"

TimeStepperTEBDnew::TimeStepperTEBDnew(const SiteSet& sites, const double J, const double tstep, const Args& args)
  : J(J), sites(sites), args(args) {
  setTstep(tstep);
}

void TimeStepperTEBDnew::initJGates(const double J){
  using Gate = BondGate<IQTensor>;

  JGates_tforwards.clear();
  JGates_tbackwards.clear();

  for(int i = 1; i < sites.N(); i += 2)
      {
      auto hterm = -J*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -J*sites.op("Adag",i)*sites.op("A",i+1);

      auto gf = Gate(sites,i,i+1,Gate::tReal,tstep,hterm);
      auto gb = Gate(sites,i,i+1,Gate::tReal,-tstep,hterm);
      JGates_tforwards.push_back(gf);
      JGates_tbackwards.push_back(gb);
      }

  int offset = 1;
  if (sites.N() % 2 == 0) { offset = 2; }

  for(int i = sites.N()-offset; i >= 1; i -= 2)
      {
      auto hterm = -J*sites.op("A",i)*sites.op("Adag",i+1);
      hterm += -J*sites.op("Adag",i)*sites.op("A",i+1);

      auto gf = Gate(sites,i,i+1,Gate::tReal,tstep,hterm);
      auto gb = Gate(sites,i,i+1,Gate::tReal,-tstep,hterm);
      JGates_tforwards.push_back(gf);
      JGates_tbackwards.push_back(gb);
      }
}


void TimeStepperTEBDnew::setTstep(const double tstep_){
  tstep = tstep_;
  initJGates(J);
}

double TimeStepperTEBDnew::getTstep(){
  return tstep;
}

void TimeStepperTEBDnew::initUGates(const double U){
  UGates.clear();

  for (int k = 1; k <= sites.N(); ++k) {
    auto s    = sites.si(k);
    auto sP   = prime(s);
    int HD    = s.nblock();

    IQTensor T(dag(s),sP);

    for (size_t i = 0; i < HD; i++) {
      T.set(s(i+1),sP(i+1), std::exp( -0.25*U*tstep*Cplx_i*i*(i-1) ) );
    }

    UGates.push_back(T);
  }
}


void TimeStepperTEBDnew::step(IQMPS& psi, const double from, const double to, const bool propagateForward){
  double U = 0.5*(from+to);

  if (propagateForward) {
    initUGates(U);
    doStep(psi, JGates_tforwards);
  }
  else {
    initUGates(-U);
    doStep(psi, JGates_tbackwards);
  }
}

void TimeStepperTEBDnew::doStep(IQMPS& psi, const GateList& JGates){
  // if N odd: "lonely" UGate at end must be applied first
  if (sites.N() % 2 != 0) { // N is odd
    auto AN = psi.A(sites.N());
    AN *= UGates.at(sites.N()-1);
    AN.mapprime(1,0,Site);
    psi.setA(sites.N(),AN);
  }

  const bool normalize = args.getBool("Normalize",true);
  bool movingFromLeft = true;
  IQTensor AA;
  auto g = JGates.begin();
  while(g != JGates.end())
      {
      auto i1 = g->i1();
      auto i2 = g->i2();
      if (movingFromLeft) {
        AA = psi.A(i1)*psi.A(i2)*UGates.at(i1-1)*UGates.at(i2-1)*prime(g->gate(),Site);

        // if N even apply lonely Ugate at right side in the end of left move
        if (i2 == sites.N() && sites.N() % 2 == 0) { // N is even
          AA *= prime(prime(UGates.at(i2-1),Site),Site);
          AA.mapprime(3,2,Site);
        }
      } else {
        AA = psi.A(i1)*psi.A(i2)*g->gate()*prime(UGates.at(i1-1),Site)*prime(UGates.at(i2-1),Site);
      }

      AA.mapprime(2,1,Site);
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
              psi.position(ni1); //does no work if position already ni1
              }
          if(ni1 < i2)
              {
              psi.svdBond(i1,AA,Fromright,args);
              psi.position(ni2); //does no work if position already ni2
              }
          if (i2 == ni1 || i1 == ni2) // odd condition || even condition
              {
              movingFromLeft = false;
              }
          }
      else
          {
          //No next gate to analyze, just restore MPS form
          psi.svdBond(i1,AA,Fromright,args);
          }
      }

  // "lonely" UGate at start must be applied last
  auto A1 = psi.A(1);
  A1 *= UGates.front();
  A1.mapprime(1,0,Site);
  psi.setA(1,A1);

  if(normalize)
      {
      psi.normalize();
      }
}
