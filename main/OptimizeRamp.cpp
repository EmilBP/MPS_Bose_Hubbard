#include "OCBoseHubbard_nlp.hpp"
#include "OptimalControl.hpp"
#include "ControlBasisFactory.hpp"
#include "SeedGenerator.hpp"
#include "IpIpoptApplication.hpp"
#include "itensor/all.h"
#include "boson.h"
#include "HamiltonianBH.hpp"
#include "TimeStepperTEBDfast.hpp"
#include "InitializeState.hpp"
#include <stdlib.h>
#include <time.h>
#include <string>

using namespace itensor;
using namespace Ipopt;


int main(){


  srand ((unsigned)time(NULL));

  double tstep  = 1e-2;
  double T      = 0.2;

  int N         = 5;
  int Npart     = 5;
  int locDim    = 5;

  double J      = 1.0;
  double U_i    = 2.0;
  double U_f    = 50;

  int M         = 8;
  double gamma  = 0;

  std::cout << "Performing optimal control of Bose-Hubbard model ... \n\n";
  std::cout << " ***** Parameters used ***** \n";
  std::cout << "Number of sites ...... " << N << "\n";
  std::cout << "Number of particles ...... " << Npart << "\n";
  std::cout << "Local Fock space dimension ...... " << locDim << "\n";
  std::cout << "Control duration ...... " << T << "\n";
  std::cout << "Time-step size ...... " << tstep << "\n";
  std::cout << "GROUP dimension ...... " << M << "\n";
  std::cout << "Gamma (regularisation factor) ...... " << gamma << "\n\n\n";

  auto sites    = Boson(N,locDim);
  auto psi_i    = InitializeState(sites,Npart,J,U_i);
  auto psi_f    = InitializeState(sites,Npart,J,U_f);

  auto H_BH     = HamiltonianBH(sites,J,tstep,0);
  auto TEBD     = TimeStepperTEBDfast(sites,J,tstep,{"Cutoff=",1E-8});
  OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, gamma);

  auto u0       = SeedGenerator::linsigmoidSeed(U_i,U_f,T/tstep+1);
  auto bControl = ControlBasisFactory::buildCBsin(u0,tstep,T,M);


  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<TNLP> mynlp = new OCBoseHubbard_nlp(OC,bControl,false);

  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-8);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  app->Options()->SetStringValue("output_file", "logfile_BH.txt");


  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return 0;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  return 0;
}
