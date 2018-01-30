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

  auto sites  = Boson(N,locDim);
  auto psi_i  = SetupSuperfluid(sites,Npart);
  auto psi_f  = SetupMottInsulator(sites,Npart);
  auto TEBD   = TimeStepperTEBD(sites,1.0,1e-2,{"Cutoff=",1E-8});

  auto ampo = AutoMPO(sites);
  for(int i = 1; i <= N; ++i) {
    ampo += 0.5,"N(N-1)",i;
  }
  auto dHdU = IQMPO(ampo);

  OptimalControl<TimeStepperTEBD> OC(psi_f,psi_i,TEBD,dHdU, 0.1 , 1e-3);

  //
  //  setup parameters
  //

  std::vector<double> control;
  for (size_t i = 0; i < 20; i++) {
    control.push_back(getRamp(i*1e-2,1));
  }

  std::cout << OC.getCost(control) << '\n';

  auto gradA = OC.getAnalyticGradient(control);
  auto gradN = OC.getNumericGradient(control);

  std::cout << "Analytic\tNumeric\tDifference\n";
  std::ofstream outFile("GradientTest.txt");
  for (size_t i = 0; i < 20; i++) {
    std::cout << gradA.at(i) << "\t" << gradN.at(i) << "\t" << gradA.at(i)-gradN.at(i) << '\n';
    outFile << gradA.at(i) << '\t';
    outFile << gradN.at(i) << '\t';
    outFile << gradA.at(i)-gradN.at(i) << '\n';
  }


  return 0;
}
