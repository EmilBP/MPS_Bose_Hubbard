#include "itensor/all.h"
#include "boson.h"
#include "TimeStepperTEBD.hpp"
#include "InitializeState.hpp"
#include "OptimalControl.hpp"
#include <fstream>

using namespace itensor;


double getRamp(double t, double T){
  return 3+11/T*t;
}

int main(){
  int N       = 5;
  int Npart   = 5;
  int locDim  = 5;

  double tstep = 4e-3;
  double J     = 1.0;

  auto sites  = Boson(N,locDim);
  auto psi_i  = SetupSuperfluid(sites,Npart);
  auto psi_f  = SetupMottInsulator(sites,Npart);
  auto TEBD   = TimeStepperTEBD(sites,J,tstep,{"Cutoff=",1E-8});

  auto ampo = AutoMPO(sites);
  auto prefac0 = Cplx_i*tstep;
  auto prefac1 = -tstep*tstep/2.0;
  auto prefac2 = -Cplx_i*tstep*tstep*tstep/6.0;
  for(int i = 1; i <= N; ++i) {
    // zeroth order
    ampo += prefac0*0.5,"N(N-1)",i;
  }
  for(int i = 1; i < N; ++i) {
    // first order
    ampo += -prefac1*2.0*J,"N",i+1,"Adag",i,"A",i+1;
    ampo += prefac1*2.0*J,"N",i,"Adag",i,"A",i+1;
    ampo += -prefac1*2.0*J,"Adag",i,"A",i+1;

    ampo += -prefac1*2.0*J,"N",i,"Adag",i+1,"A",i;
    ampo += prefac1*2.0*J,"N",i+1,"Adag",i+1,"A",i;
    ampo += -prefac1*2.0*J,"Adag",i+1,"A",i;
  }
    // // second order
    // // [H_J , [H,H_U]]
    // ampo += prefac2*6.0*J*J,"N",i+1,"N",i+1;
    // ampo += -prefac2*4.0*J*J,"N",i+1,"N",i;
    // ampo += -prefac2*2.0*J*J,"N",i+1,"N",i+2;
    // ampo += prefac2*2.0*J*J,"N",i+1,"Adag",i+1,"A",i+3;
    // ampo += prefac2*2.0*J*J,"N",i+1,"Adag",i+3,"A",i+1;
    // ampo += -prefac2*4.0*J*J,"N",i+1,"Adag",i+2,"A",i;
    // ampo += -prefac2*4.0*J*J,"N",i+1,"Adag",i,"A",i+2;
    // ampo += prefac2*2.0*J*J,"N",i+1,"Adag",i-1,"A",i+1;
    // ampo += prefac2*2.0*J*J,"N",i+1,"Adag",i+1,"A",i-1;
    //
    // ampo += prefac2*4.0*J*J,"Adag",i,"A",i+1 , "Adag",i,"A",i+1;
    // ampo += -prefac2*4.0*J*J,"Adag",i,"A",i+1 , "Adag",i+1,"A",i;
    // ampo += prefac2*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i+2,"A",i+1;
    // ampo += -prefac2*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i+1,"A",i+2;
    // ampo += prefac2*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i,"A",i-1;
    // ampo += -prefac2*2.0*J*J,"Adag",i,"A",i+1 , "Adag",i-1,"A",i;
    //
    // ampo += prefac2*4.0*J*J,"Adag",i+1,"A",i , "Adag",i+1,"A",i;
    // ampo += -prefac2*4.0*J*J,"Adag",i+1,"A",i , "Adag",i,"A",i+1;
    // ampo += prefac2*2.0*J*J,"Adag",i+1,"A",i , "Adag",i+1,"A",i+2;
    // ampo += -prefac2*2.0*J*J,"Adag",i+1,"A",i , "Adag",i+2,"A",i+1;
    // ampo += prefac2*2.0*J*J,"Adag",i+1,"A",i , "Adag",i-1,"A",i;
    // ampo += -prefac2*2.0*J*J,"Adag",i+1,"A",i , "Adag",i,"A",i-1;
    //
    // // [H_U , [H,H_U]]
    // ampo += -prefac2*2.0*J*U,"N",i+1,"N",i,"Adag",i,"A",i+1;
    // ampo += prefac2*2.0*J*U,"N",i+1,"N",i+1,"Adag",i,"A",i+1;
    // ampo += prefac2*2.0*J*U,"N",i+1,"Adag",i,"A",i+1;
    //
    // ampo += prefac2*2.0*J*U,"N",i,"N",i,"Adag",i,"A",i+1;
    // ampo += -prefac2*2.0*J*U,"N",i,"N",i+1,"Adag",i,"A",i+1;
    // ampo += -prefac2*2.0*J*U,"N",i,"Adag",i,"A",i+1;
    //
    // ampo += -prefac2*2.0*J*U,"N",i,"Adag",i,"A",i+1;
    // ampo += prefac2*2.0*J*U,"N",i+1,"Adag",i,"A",i+1;
    // ampo += prefac2*2.0*J*U,"Adag",i,"A",i+1;
    //
    //
    // ampo += -prefac2*2.0*J*U,"N",i,"N",i+1,"Adag",i+1,"A",i;
    // ampo += prefac2*2.0*J*U,"N",i,"N",i,"Adag",i+1,"A",i;
    // ampo += prefac2*2.0*J*U,"N",i,"Adag",i+1,"A",i;
    //
    // ampo += prefac2*2.0*J*U,"N",i+1,"N",i+1,"Adag",i+1,"A",i;
    // ampo += -prefac2*2.0*J*U,"N",i+1,"N",i,"Adag",i+1,"A",i;
    // ampo += -prefac2*2.0*J*U,"N",i+1,"Adag",i+1,"A",i;
    //
    // ampo += -prefac2*2.0*J*U,"N",i+1,"Adag",i+1,"A",i;
    // ampo += prefac2*2.0*J*U,"N",i,"Adag",i+1,"A",i;
    // ampo += prefac2*2.0*J*U,"Adag",i+1,"A",i;

  auto dHdU = IQMPO(ampo);

  OptimalControl<TimeStepperTEBD> OC(psi_f,psi_i,TEBD,dHdU, 0);

  //
  //  setup parameters
  //

  double controlLength = 2.5*1e2;

  std::vector<double> control;
  for (size_t i = 0; i < controlLength; i++) {
    control.push_back(getRamp(i*tstep,controlLength*tstep));
  }

  // auto fidtest = OC.checkFidelity(control);
  // for (auto& fid : fidtest){
  //   std::cout << fid << '\n';
  // }


  auto gradA = OC.getAnalyticGradient(control);
  // auto gradN = OC.getNumericGradient(control);

  std::cout << "Analytic\tNumeric\tDifference\n";
  // std::ofstream outFile1("AnalyticGradient.txt");
  // std::ofstream outFile2("NumericGradient_tstep0.001_controlL1000.txt");
  std::fstream myfile("NumericGradient_tstep0.001_controlL1000.txt", std::ios_base::in);

  double gN;
  size_t count = 0;
  while (myfile >> gN){
    std::cout << gradA.second.at(count) << "\t";
    std::cout << gN << "\t";
    std::cout << gradA.second.at(count)-gN << "\n";
    count++;
  }

  return 0;
}
