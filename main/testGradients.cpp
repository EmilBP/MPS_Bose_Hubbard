#include "itensor/all.h"
#include "boson.h"
#include "OptimalControl.hpp"
#include "InitializeState.hpp"
#include "BoseHubbardMPO.h"
#include <math.h>

using namespace itensor;


double getRamp(double t){
  return 3+14/0.1*t;
}

int main(){
  int N       = 5;
  int Npart   = 5;
  int locDim  = 5;

  auto sites  = Boson(N,locDim);
  auto psi_i  = SetupSuperfluid(sites,Npart);
  auto psi_f  = SetupMottInsulator(sites,Npart);
  auto BHMPO  = BoseHubbardMPO(sites);
  auto gamma  = 0;

  auto opt    = OptimalControl(psi_f, psi_i, BHMPO, gamma);
  
  //
  //  setup parameters
  //

  auto args   = Args("Cutoff=",1E-11,"Maxm=",400);
  double eps  = 1e-6;
  double dt   = 1e-2;
  double T    = 0.1;
  double temp = 0;

  std::vector<double> times;
  std::vector<double> ramp;

  times.emplace_back(temp);
  ramp.emplace_back(getRamp(temp));
  while (T-temp > 1e-6) {
    temp += dt;
    times.emplace_back(temp);
    ramp.emplace_back(getRamp(temp));
  }

  for (auto& u : ramp){
    std::cout << "V0: " << u << " J: " << BHMPO.getJ(u) << ". U: " << BHMPO.getU(u) << std::endl;
  }

  auto anal   = opt.getAnalyticGradientTest(ramp,dt,args);
  auto nums   = opt.getNumericGradient(ramp,eps,dt,args);

  std::cout << anal.size() << std::endl;
  std::cout << nums.size() << std::endl;
  std::cout << ramp.size() << std::endl;
  std::cout << "Analytical \t\t Numerical \t\t Difference" << std::endl;
  for (size_t i = 0; i < nums.size(); i++) {
    std::cout << anal[i] << "\t\t" << nums[i] << "\t\t" << anal[i]-nums[i] << "\n";
  }

  return 0;
}
