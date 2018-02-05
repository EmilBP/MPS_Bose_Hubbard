#include "HamiltonianBH.hpp"

HamiltonianBH::HamiltonianBH(const SiteSet& sites, const double J)
  : sites(sites), J(J), N(sites.N()) {

}

IQMPO HamiltonianBH::dHdU(const double& control_n, const double& tstep){
  auto ampo     = AutoMPO(sites);
  auto prefac0  = Cplx_i*tstep;
  auto prefac1  = -tstep*tstep/2.0;
  auto prefac2  = -Cplx_i*tstep*tstep*tstep/6.0;

  for(int i = 1; i <= N; ++i) {
    ampo += prefac0*0.5,"N(N-1)",i;

    ampo += prefac2*6.0*J*J,"N",i,"N",i;
  }

  for(int i = 1; i <= N-1; ++i) {
    ampo += -prefac1*2.0*J,"N",i+1,"Adag",i,"A",i+1;
    ampo += prefac1*2.0*J,"N",i,"Adag",i,"A",i+1;
    ampo += -prefac1*2.0*J,"Adag",i,"A",i+1;
    ampo += -prefac1*2.0*J,"N",i,"Adag",i+1,"A",i;
    ampo += prefac1*2.0*J,"N",i+1,"Adag",i+1,"A",i;
    ampo += -prefac1*2.0*J,"Adag",i+1,"A",i;

    // [H_J , [H,H_U]]
    ampo += -prefac2*4.0*J*J,"N",i+1,"N",i;
    ampo += -prefac2*2.0*J*J,"N",i,"N",i+1;
    ampo += prefac2*4.0*J*J,"Adag",i,"A",i+1 , "Adag",i,"A",i+1;
    ampo += -prefac2*4.0*J*J,"Adag",i,"A",i+1 , "Adag",i+1,"A",i;
    ampo += prefac2*4.0*J*J,"Adag",i+1,"A",i , "Adag",i+1,"A",i;
    ampo += -prefac2*4.0*J*J,"Adag",i+1,"A",i , "Adag",i,"A",i+1;

    // [H_U , [H,H_U]]
    ampo += -prefac2*2.0*J*control_n,"N",i+1,"N",i,"Adag",i,"A",i+1;
    ampo += prefac2*2.0*J*control_n,"N",i+1,"N",i+1,"Adag",i,"A",i+1;
    ampo += prefac2*2.0*J*control_n,"N",i+1,"Adag",i,"A",i+1;
    ampo += prefac2*2.0*J*control_n,"N",i,"N",i,"Adag",i,"A",i+1;
    ampo += -prefac2*2.0*J*control_n,"N",i,"N",i+1,"Adag",i,"A",i+1;
    ampo += -prefac2*2.0*J*control_n,"N",i,"Adag",i,"A",i+1;
    ampo += -prefac2*2.0*J*control_n,"N",i,"Adag",i,"A",i+1;
    ampo += prefac2*2.0*J*control_n,"N",i+1,"Adag",i,"A",i+1;
    ampo += prefac2*2.0*J*control_n,"Adag",i,"A",i+1;
    ampo += -prefac2*2.0*J*control_n,"N",i,"N",i+1,"Adag",i+1,"A",i;
    ampo += prefac2*2.0*J*control_n,"N",i,"N",i,"Adag",i+1,"A",i;
    ampo += prefac2*2.0*J*control_n,"N",i,"Adag",i+1,"A",i;
    ampo += prefac2*2.0*J*control_n,"N",i+1,"N",i+1,"Adag",i+1,"A",i;
    ampo += -prefac2*2.0*J*control_n,"N",i+1,"N",i,"Adag",i+1,"A",i;
    ampo += -prefac2*2.0*J*control_n,"N",i+1,"Adag",i+1,"A",i;
    ampo += -prefac2*2.0*J*control_n,"N",i+1,"Adag",i+1,"A",i;
    ampo += prefac2*2.0*J*control_n,"N",i,"Adag",i+1,"A",i;
    ampo += prefac2*2.0*J*control_n,"Adag",i+1,"A",i;
  }

  for(int i = 1; i <= N-2; ++i) {
    // [H_J , [H,H_U]]
    ampo += prefac2*2.0*J*J,"N",i,"Adag",i,"A",i+2;
    ampo += prefac2*2.0*J*J,"N",i,"Adag",i+2,"A",i;
    ampo += -prefac2*4.0*J*J,"N",i+1,"Adag",i+2,"A",i;
    ampo += -prefac2*4.0*J*J,"N",i+1,"Adag",i,"A",i+2;
    ampo += prefac2*2.0*J*J,"N",i+2,"Adag",i,"A",i+2;
    ampo += prefac2*2.0*J*J,"N",i+2,"Adag",i+2,"A",i;

    ampo += prefac2*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i+2,"A",i+1;
    ampo += -prefac2*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i+1,"A",i+2;
    ampo += prefac2*2.0*J*J,"Adag",i+1,"A",i+2 , "Adag",i+1,"A",i;
    ampo += -prefac2*2.0*J*J,"Adag",i+1,"A",i+2 , "Adag",i,"A",i+1;

    ampo += prefac2*2.0*J*J,"Adag",i+1,"A",i , "Adag",i+1,"A",i+2;
    ampo += -prefac2*2.0*J*J,"Adag",i+1,"A",i , "Adag",i+2,"A",i+1;
    ampo += prefac2*2.0*J*J,"Adag",i+2,"A",i+1 , "Adag",i,"A",i+1;
    ampo += -prefac2*2.0*J*J,"Adag",i+2,"A",i+1 , "Adag",i+1,"A",i;
  }

  return IQMPO(ampo);
}
