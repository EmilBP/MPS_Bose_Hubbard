#include "itensor/all.h"
#include "boson.h"
#include "correlations.h"
#include "OptimalControl.hpp"
#include "BoseHubbardMPO.h"
#include "InitializeState.hpp"
#include <vector>
#include <fstream>
#include <ctime>

using namespace itensor;
using std::vector;

int main()
  {
  int N       = 5;
  int Npart   = 5;
  int locDim  = 5;

  auto sites  = Boson(N,locDim);
  auto psi_i  = SetupSuperfluid(sites,Npart);
  auto psi_f  = SetupMottInsulator(sites,Npart);
  auto BHMPO  = BoseHubbardMPO(sites,0,0,0);
  auto gamma  = 0;

  auto opt    = OptimalControl(psi_f, psi_i, BHMPO, gamma);

  //
  //  setup parameters
  //

  auto args       = Args("Cutoff=",1E-9,"Maxm=",50);
  double dt       = 1e-2;
  auto T          = 1;
  size_t Nsteps   = T/dt;
  size_t maxeval  = 1e2;
  vector<double> control(Nsteps, 0.0);

  vector<double> costs = opt.OptimizeControl(control, dt, maxeval, args);

  return 0;
  }
