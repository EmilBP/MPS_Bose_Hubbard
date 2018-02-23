#include "HamiltonianBH.hpp"

HamiltonianBH::HamiltonianBH(const SiteSet& sites, double J, double tstep, size_t expansionOrder)
  : sites(sites), J(J), tstep(tstep), N(sites.N()), ampo(AutoMPO(sites)) {

  setTstep(tstep);
  setExpansionOrder(expansionOrder);
}

void HamiltonianBH::setTstep(const double tstep_){
  tstep = tstep_;
  prefactors.clear();

  prefactors.push_back(tstep);
  prefactors.push_back(Cplx_i*tstep*tstep/2.0);
  prefactors.push_back(-tstep*tstep*tstep/6.0);
}

void HamiltonianBH::setExpansionOrder(size_t expansionOrder_){
  expansionOrder = expansionOrder_;

  ampo = AutoMPO(sites);
  for(int i = 1; i <= N; ++i) {
    ampo += prefactors.at(0)*0.5,"N(N-1)",i;

    if (expansionOrder >= 2){
      ampo += prefactors.at(2)*3.0*J*J,"N",i,"N",i;
    }
  }

  for(int i = 1; i <= N-1; ++i) {
    if (expansionOrder >= 1) {
      ampo += -prefactors.at(1)*J,"N",i+1,"Adag",i,"A",i+1;
      ampo += prefactors.at(1)*J,"N",i,"Adag",i,"A",i+1;
      ampo += -prefactors.at(1)*J,"Adag",i,"A",i+1;
      ampo += -prefactors.at(1)*J,"N",i,"Adag",i+1,"A",i;
      ampo += prefactors.at(1)*J,"N",i+1,"Adag",i+1,"A",i;
      ampo += -prefactors.at(1)*J,"Adag",i+1,"A",i;
      // ampo += -prefactors.at(1)*0.5*J,"N",i+1,"Adag",i,"A",i+1;
      // ampo += prefactors.at(1)*0.5*J,"N",i,"Adag",i,"A",i+1;
      // ampo += -prefactors.at(1)*0.5*J,"Adag",i,"A",i+1,"N",i+1;
      // ampo += prefactors.at(1)*0.5*J,"Adag",i,"A",i+1,"N",i;
      // ampo += -prefactors.at(1)*0.5*J,"N",i,"Adag",i+1,"A",i;
      // ampo += prefactors.at(1)*0.5*J,"N",i+1,"Adag",i+1,"A",i;
      // ampo += -prefactors.at(1)*0.5*J,"Adag",i+1,"A",i,"N",i;
      // ampo += prefactors.at(1)*0.5*J,"Adag",i+1,"A",i,"N",i+1;
    }
    if (expansionOrder >= 2) {
      // [H_J , [H,H_U]]
      ampo += -prefactors.at(2)*2.0*J*J,"N",i+1,"N",i;
      ampo += -prefactors.at(2)*2.0*J*J,"N",i,"N",i+1; // <-- check this
      ampo += prefactors.at(2)*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i,"A",i+1;
      ampo += -prefactors.at(2)*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i+1,"A",i;
      ampo += prefactors.at(2)*2.0*J*J,"Adag",i+1,"A",i , "Adag",i+1,"A",i;
      ampo += -prefactors.at(2)*2.0*J*J,"Adag",i+1,"A",i , "Adag",i,"A",i+1;
    }
  }

  if (expansionOrder >= 2) {
    for(int i = 1; i <= N-2; ++i) {
      // [H_J , [H,H_U]]
      ampo += prefactors.at(2)*J*J,"N",i,"Adag",i,"A",i+2;
      ampo += prefactors.at(2)*J*J,"N",i,"Adag",i+2,"A",i;
      ampo += -prefactors.at(2)*2.0*J*J,"N",i+1,"Adag",i+2,"A",i;
      ampo += -prefactors.at(2)*2.0*J*J,"N",i+1,"Adag",i,"A",i+2;
      ampo += prefactors.at(2)*J*J,"N",i+2,"Adag",i,"A",i+2;
      ampo += prefactors.at(2)*J*J,"N",i+2,"Adag",i+2,"A",i;

      ampo += prefactors.at(2)*J*J,"Adag",i,"A",i+1 , "Adag",i+2,"A",i+1;
      ampo += -prefactors.at(2)*J*J,"Adag",i,"A",i+1 , "Adag",i+1,"A",i+2;
      ampo += prefactors.at(2)*J*J,"Adag",i+1,"A",i+2 , "Adag",i+1,"A",i;
      ampo += -prefactors.at(2)*J*J,"Adag",i+1,"A",i+2 , "Adag",i,"A",i+1;

      ampo += prefactors.at(2)*J*J,"Adag",i+1,"A",i , "Adag",i+1,"A",i+2;
      ampo += -prefactors.at(2)*J*J,"Adag",i+1,"A",i , "Adag",i+2,"A",i+1;
      ampo += prefactors.at(2)*J*J,"Adag",i+2,"A",i+1 , "Adag",i,"A",i+1;
      ampo += -prefactors.at(2)*J*J,"Adag",i+2,"A",i+1 , "Adag",i+1,"A",i;
    }
  }

  if (expansionOrder < 2) {
    dHdUconst = IQMPO(ampo);
  }
}


IQMPO HamiltonianBH::dHdU(const double& control_n){
  if (expansionOrder < 2) {
    return dHdUconst;
  }
  else {
    for(int i = 1; i <= N-1; ++i) {
      // [H_U , [H,H_U]]
      ampo += -prefactors.at(2)*J*control_n,"N",i+1,"N",i,"Adag",i,"A",i+1;
      ampo += prefactors.at(2)*J*control_n,"N",i+1,"N",i+1,"Adag",i,"A",i+1;
      ampo += prefactors.at(2)*J*control_n,"N",i+1,"Adag",i,"A",i+1;
      ampo += prefactors.at(2)*J*control_n,"N",i,"N",i,"Adag",i,"A",i+1;
      ampo += -prefactors.at(2)*J*control_n,"N",i,"N",i+1,"Adag",i,"A",i+1;
      ampo += -prefactors.at(2)*J*control_n,"N",i,"Adag",i,"A",i+1;
      ampo += -prefactors.at(2)*J*control_n,"N",i,"Adag",i,"A",i+1;
      ampo += prefactors.at(2)*J*control_n,"N",i+1,"Adag",i,"A",i+1;
      ampo += prefactors.at(2)*J*control_n,"Adag",i,"A",i+1;
      ampo += -prefactors.at(2)*J*control_n,"N",i,"N",i+1,"Adag",i+1,"A",i;
      ampo += prefactors.at(2)*J*control_n,"N",i,"N",i,"Adag",i+1,"A",i;
      ampo += prefactors.at(2)*J*control_n,"N",i,"Adag",i+1,"A",i;
      ampo += prefactors.at(2)*J*control_n,"N",i+1,"N",i+1,"Adag",i+1,"A",i;
      ampo += -prefactors.at(2)*J*control_n,"N",i+1,"N",i,"Adag",i+1,"A",i;
      ampo += -prefactors.at(2)*J*control_n,"N",i+1,"Adag",i+1,"A",i;
      ampo += -prefactors.at(2)*J*control_n,"N",i+1,"Adag",i+1,"A",i;
      ampo += prefactors.at(2)*J*control_n,"N",i,"Adag",i+1,"A",i;
      ampo += prefactors.at(2)*J*control_n,"Adag",i+1,"A",i;
    }

    return IQMPO(ampo);
  }
}
