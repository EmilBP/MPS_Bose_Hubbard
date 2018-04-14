#include "OptimalControl.hpp"
#include "SeedGenerator.hpp"
#include "itensor/all.h"
#include "boson.h"
#include "HamiltonianBH.hpp"
#include "TimeStepperTEBDfast.hpp"
#include "InitializeState.hpp"
#include <stdlib.h>
#include <time.h>
#include <string>

using namespace itensor;

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

  double gamma  = 0;
  int dHorder   = 0;
  int seed      = 1;

  if(argc > 2) seed = std::stoi(argv[2]);
  else printfln("Default seed used");

  srand ((unsigned) seed*time(NULL));


  std::cout << "Running adiabatic ramp for Bose-Hubbard model ... \n\n";
  std::cout << " ******* Parameters used ******* \n";
  std::cout << "Number of sites ................ " << N << "\n";
  std::cout << "Number of particles ............ " << Npart << "\n";
  std::cout << "Local Fock space dimension ..... " << locDim << "\n";
  std::cout << "Control duration ............... " << T << "\n";
  std::cout << "Time-step size ................. " << tstep << "\n";
  std::cout << "Gamma (regularisation) ......... " << gamma << "\n";
  std::cout << "Seed  .......................... " << seed << "\n\n\n";


  auto sites    = Boson(N,locDim);
  auto psi_i    = InitializeState(sites,Npart,J,U_i);
  auto psi_f    = InitializeState(sites,Npart,J,U_f);

  auto H_BH     = HamiltonianBH(sites,J,tstep,dHorder);
  auto TEBD     = TimeStepperTEBDfast(sites,J,tstep,{"Cutoff=",1E-8});
  OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH,gamma);

  auto times  = SeedGenerator::generateRange(0,tstep,T);
  auto u0     = SeedGenerator::linsigmoidSeed(U_i,U_f,T/tstep+1);
  auto F      = OC.getFidelityForAllT(u0);


  std::string filename = "AdiabaticF.txt";
  std::ofstream myfile (filename);
  if (myfile.is_open())
  {
    for (size_t i = 0; i < times.size(); i++) {
      myfile << times.at(i) << "\t";
      myfile << u0.at(i) << "\t";
      myfile << F.at(i) << "\n";
    }
    myfile.close();
  }
  else std::cout << "Unable to open file\n";


  return 0;
}
