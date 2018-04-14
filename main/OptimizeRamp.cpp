#include "OCBoseHubbard_nlp.hpp"
#include "OptimalControl.hpp"
#include "ControlBasisFactory.hpp"
#include "SeedGenerator.hpp"
#include "IpIpoptApplication.hpp"
#include "itensor/all.h"
#include "boson.h"
#include "correlations.hpp"
#include "HamiltonianBH.hpp"
#include "TimeStepperTEBDfast.hpp"
#include "InitializeState.hpp"
#include <stdlib.h>
#include <time.h>
#include <string>

using namespace itensor;
using namespace Ipopt;


int main(int argc, char* argv[]){

  if(argc < 2) {
    printfln("Usage: %s InputFile_BHcontrol",argv[0]);
    return 0;
  }

  auto input    = InputGroup(argv[1],"input");

  double tstep  = input.getReal("tstep",1e-2);
  double T      = input.getReal("T");

  int N         = input.getInt("N");
  int Npart     = input.getInt("Npart");
  int locDim    = input.getInt("d");

  double J      = 1.0;
  double U_i    = 2.0;
  double U_f    = 50;

  int M         = input.getInt("M");
  double gamma  = input.getReal("gamma",0);
  int dHorder   = input.getInt("dHOrder",0);
  bool cache    = input.getYesNo("cacheProgress",false);
  int seed      = 1;

  if(argc > 2) seed = std::stoi(argv[2]);
  else printfln("Default seed used");

  srand ((unsigned) seed*time(NULL));


  std::cout << "Performing optimal control of Bose-Hubbard model ... \n\n";
  std::cout << " ******* Parameters used ******* \n";
  std::cout << "Number of sites ................ " << N << "\n";
  std::cout << "Number of particles ............ " << Npart << "\n";
  std::cout << "Local Fock space dimension ..... " << locDim << "\n";
  std::cout << "Control duration ............... " << T << "\n";
  std::cout << "Time-step size ................. " << tstep << "\n";
  std::cout << "GROUP dimension ................ " << M << "\n";
  std::cout << "Gamma (regularisation) ......... " << gamma << "\n";
  std::cout << "Seed  .......................... " << seed << "\n\n\n";


  auto sites    = Boson(N,locDim);
  auto psi_i    = InitializeState(sites,Npart,J,U_i);
  auto psi_f    = InitializeState(sites,Npart,J,U_f);

  auto H_BH     = HamiltonianBH(sites,J,tstep,dHorder);
  auto TEBD     = TimeStepperTEBDfast(sites,J,tstep,{"Cutoff=",1E-8});
  auto times    = SeedGenerator::generateRange(0,tstep,T);
  OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH,gamma);

  auto u0       = SeedGenerator::linsigmoidSeed(U_i,U_f,T/tstep+1);
  auto bControl = ControlBasisFactory::buildCBsin(u0,tstep,T,M);

  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<TNLP> mynlp = new OCBoseHubbard_nlp(OC,bControl,times,cache);

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
  // app->Options()->SetStringValue("output_file", "logfile_BH.txt");


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

  // Extract psi for each t, evaluate expectation value
  // of number operator, and save to file.
  auto psi_t = OC.getPsit();
  std::string filename = "ExpectationN.txt";
  std::ofstream myfile (filename);
  if (myfile.is_open())
  {
    size_t ind = 0;
    for (auto& psi : psi_t){
      myfile << times.at(ind++) << "\t";
      auto expn = expectationValues(sites,psi,"N");
      for (auto& val : expn){
        myfile << val.real() << "\t";
      }
      myfile << "\n";
    }
    myfile.close();
  }
  else std::cout << "Unable to open file\n";


  return 0;
}
