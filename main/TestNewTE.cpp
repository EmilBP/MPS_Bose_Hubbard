#include "itensor/all.h"
#include "boson.h"
#include "TimeStepperTEBD.hpp"
#include "HamiltonianBH.hpp"
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
  auto H_BH   = HamiltonianBH(sites,J);

  OptimalControl<TimeStepperTEBD,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);

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
