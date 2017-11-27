#include "itensor/all.h"
#include "boson.h"
#include "OptimalControl.hpp"
#include "InitializeState.hpp"
#include "BoseHubbardMPO.h"
#include <math.h>

using namespace itensor;


double getRamp(double t){
  // double T = 0.1;
  // return 8.0*pow(sin(2*M_PI*t/T),2.0);
  return exp(27.08*t)-1;
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

  auto args   = Args("Cutoff=",1E-9,"Maxm=",50);
  double eps  = 1e-4;
  double dt   = 1e-3;
  double T    = 0.1;
  double temp = 0;

  std::vector<double> times;
  std::vector<double> ramp;

  times.emplace_back(temp);
  ramp.emplace_back(getRamp(temp));
  while (temp < T) {
    temp += dt;
    times.emplace_back(temp);
    ramp.emplace_back(getRamp(temp));
  }

  auto anal   = opt.getAnalyticGradient(ramp,dt);
  std::cout << "hej" << std::endl; 
  auto nums   = opt.getNumericGradient(ramp,eps,dt,args);

  for (size_t i = 0; i < nums.size(); i++) {
    std::cout << "Analytical gradient: " << anal[i] << ". Numerical gradient: " << nums[i] << ". Difference: " << anal[i]-nums[i] << std::endl;
  }

  return 0;
}
