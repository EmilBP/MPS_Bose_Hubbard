#include "OCdummy_nlp.hpp"
#include "OptimalControlDummy.hpp"
#include "OptimalControl.hpp"
#include "ControlBasisFactory.hpp"
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

using matrix = std::vector< std::vector<double> >;

std::vector<double> generateRange(double a, double b, double c) { //equiv to a:b:c
    std::vector<double> array;
    while(a <= c + 1e-7) {
        array.push_back(a);
        a += b;         // could recode to better handle rounding errors
    }
    return array;
}

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> array;
    double step = (b-a) / (n-1);

    while(a <= b + 1e-7) {
        array.push_back(a);
        a += step;           // could recode to better handle rounding errors
    }
    return array;
}

double randomDouble(double min, double max){
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}

void printData(const matrix& data){
  for (auto& row : data){
    for (auto& val : row){
      std::cout << val << "\t";
    }
    std::cout << "\n";
  }
}

std::vector<double> randomVec(double min, double max, size_t N){
  std::vector<double> v;
  for (size_t i = 0; i < N; i++) {
    v.push_back(randomDouble(min,max));
  }
  return v;
}

void saveData(const matrix& data, const std::string filename){
  std::ofstream myfile (filename);
  if (myfile.is_open())
  {
    for (auto& row : data){
      for (auto& val : row){
        myfile << val << "\t";
      }
      myfile << "\n";
    }
    myfile.close();
  }
  else std::cout << "Unable to open file\n";
}


matrix matchGradients(std::vector<double> weights, double tstep, double cstart, double cend, double T) {
  matrix result;

  auto OCD      = OptimalControlDummy(weights,tstep);
  auto times    = generateRange(0,tstep,T);
  auto control  = linspace(cstart,cend,times.size());
  auto Agrad    = OCD.getAnalyticGradient(control);
  auto Ngrad    = OCD.getNumericGradient(control);

  for (size_t i = 0; i < Agrad.second.size(); i++) {
    std::vector<double> tmp;
    tmp.push_back(Agrad.second.at(i));
    tmp.push_back(Ngrad.second.at(i));
    result.push_back(tmp);
  }

  return result;
}

matrix matchControlGradients(std::vector<double> weights, double tstep, double cstart, double cend, double T) {
  matrix result;

  auto OCD      = OptimalControlDummy(weights,tstep);
  auto c        = linspace(cstart,cend,10);
  auto u0       = linspace(0,10,T/tstep+1);
  auto control  = ControlBasisFactory::buildCBsin(u0,tstep,T,c.size());

  auto Agrad    = OCD.getAnalyticGradient(control);
  auto Ngrad    = OCD.getNumericGradient(control);

  for (size_t i = 0; i < Agrad.second.size(); i++) {
    std::vector<double> tmp;
    tmp.push_back(Agrad.second.at(i));
    tmp.push_back(Ngrad.second.at(i));
    result.push_back(tmp);
  }

  return result;
}

matrix matchControlGradientsBH( double tstep, double T, size_t M) {
  matrix result;
  int N         = 5;
  int Npart     = 5;
  int locDim    = 5;
  double J      = 1.0;

  auto sites    = Boson(N,locDim);
  auto psi_i    = InitializeState(sites,Npart,J,5.0);
  auto psi_f    = InitializeState(sites,Npart,J,20.0);

  auto H_BH     = HamiltonianBH(sites,J,tstep,0);
  auto TEBD     = TimeStepperTEBDfast(sites,J,tstep,{"Cutoff=",1E-8});
  OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);

  auto c        = randomVec(-2,2,M);
  auto u0       = linspace(5.0,20.0,T/tstep+1);
  auto control  = ControlBasisFactory::buildCBsin(u0,tstep,T,c.size());
  control.setCArray(c);

  auto Ngrad    = OC.getNumericGradient(control);
  auto Agrad0   = OC.getAnalyticGradient(control);
  H_BH.setExpansionOrder(1);
  auto Agrad1   = OC.getAnalyticGradient(control);

  for (size_t i = 0; i < Ngrad.second.size(); i++) {
    std::vector<double> tmp;
    tmp.push_back(Ngrad.second.at(i));
    tmp.push_back(Agrad0.second.at(i));
    tmp.push_back(Agrad1.second.at(i));

    result.push_back(tmp);
  }

  return result;
}

void testCBsinusParametrization(double tstep, double T, size_t M){
  auto u0       = linspace(0,10,T/tstep+1);
  auto bControl = ControlBasisFactory::buildCBsin(u0,tstep,T,M);

  std::vector<double> c;
  for (size_t i = 0; i < M; i++) {
    c.push_back(randomDouble(-5,5));
  }

  bControl.setCArray(c);

  std::vector<double> S,u;
  std::vector<std::vector<double> > f;
  bControl.exportParameters(u,u0,S,c,f);

  matrix USdata;
  USdata.push_back(u);
  USdata.push_back(u0);
  USdata.push_back(S);

  f.push_back(c);
  std::string name1 = "CBsinData_US_M" + std::to_string(M) + ".txt";
  std::string name2 = "CBsinData_FC_M" + std::to_string(M) + ".txt";
  saveData(USdata,name1);
  saveData(f,name2);
}


int main(){

  double tstep  = 1e-2;
  double T      = 1;
  double cstart = 2;
  double cend   = 7;

  // std::vector<int> Ms = {1, 3, 5, 7};
  // srand ((unsigned)time(NULL));
  // for (auto& M : Ms){
  //   testCBsinusParametrization(tstep,T,M);
  // }

  std::vector<size_t> Ms = {5,10,15,20,25};
  srand ((unsigned)time(NULL));
  for (auto& M : Ms){
    auto data = matchControlGradientsBH(tstep,T,M);
    std::string name = "ControlGradients_tstep001_M" + std::to_string(M) + ".txt";
    saveData(data,name);
  }


  // std::vector<double> weights = {5.5 , 1.2 , 6.3 , 0.3};
  // // auto data = matchControlGradients(weights,tstep,cstart,cend,T);
  // // printData(data);
  //
  // auto OCD      = OptimalControlDummy(weights,tstep);
  // auto u0       = linspace(0,10,T/tstep+1);
  // auto bControl = ControlBasisFactory::buildCBsin(u0,tstep,T,10);
  //
  //
  // // Create a new instance of your nlp
  // //  (use a SmartPtr, not raw)
  // SmartPtr<TNLP> mynlp = new OCdummy_nlp(OCD,bControl);
  //
  // // Create a new instance of IpoptApplication
  // //  (use a SmartPtr, not raw)
  // // We are using the factory, since this allows us to compile this
  // // example with an Ipopt Windows DLL
  // SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  //
  // // Change some options
  // // Note: The following choices are only examples, they might not be
  // //       suitable for your optimization problem.
  // app->Options()->SetNumericValue("tol", 1e-9);
  // app->Options()->SetStringValue("mu_strategy", "adaptive");
  // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  //
  // // Intialize the IpoptApplication and process the options
  // ApplicationReturnStatus status;
  // status = app->Initialize();
  // if (status != Solve_Succeeded) {
  //   printf("\n\n*** Error during initialization!\n");
  //   return (int) status;
  // }
  //
  // // Ask Ipopt to solve the problem
  // status = app->OptimizeTNLP(mynlp);
  //
  // if (status == Solve_Succeeded) {
  //   printf("\n\n*** The problem solved!\n");
  // }
  // else {
  //   printf("\n\n*** The problem FAILED!\n");
  // }
  //
  // // As the SmartPtrs go out of scope, the reference count
  // // will be decremented and the objects will automatically
  // // be deleted.


  return 0;
}
