#include "itensor/all.h"
#include "boson.h"
#include "HamiltonianBH.hpp"
#include "TimeStepperTEBD.hpp"
#include "InitializeState.hpp"
#include "OptimalControl.hpp"
#include <fstream>
#include <string>

using namespace itensor;

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

double calculateVariance(const std::vector<double>& V){
  double sum = 0, var = 0;

  for (auto& val : V){
    sum += val;
  }
  double avg = sum/V.size();
  for (auto& val : V){
    var += (val-avg)*(val-avg);
  }
  return var/V.size();
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

matrix
testCostPlusFidelity( SiteSet& sites,
                      IQMPS& psi_i,
                      IQMPS& psi_f,
                      std::vector<double>& tsteps,
                      double cstart,
                      double cend,
                      double T,
                      double J)
{

  matrix result;
  size_t count  = 1;
  auto H_BH     = HamiltonianBH(sites,J,0);

  for (auto& dt: tsteps){
    auto times    = generateRange(0,dt,T);
    auto control  = linspace(cstart,cend,times.size());

    std::cout << "Calculating using control of length " << control.size() << " w. final value " << control.back() << '\n';

    auto TEBD     = TimeStepperTEBD(sites,J,dt,{"Cutoff=",1E-9});
    OptimalControl<TimeStepperTEBD,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);
    auto CF       = OC.checkCostPlusFidelity(control);
    double var    = calculateVariance(CF.second);

    std::vector<double> tmp;
    tmp.push_back(dt);
    tmp.push_back(CF.first);
    tmp.push_back(var);

    result.push_back(tmp);
    std::cout << "Done " << count++ << " out of " << tsteps.size() << '\n';
  }

  return result;
}

matrix
matchGradients(       SiteSet& sites,
                      IQMPS& psi_i,
                      IQMPS& psi_f,
                      double tstep,
                      double cstart,
                      double cend,
                      double T,
                      double J,
                      size_t order)
{

  matrix result;
  size_t count  = 1;
  auto H_BH     = HamiltonianBH(sites,J,tstep,order);
  auto TEBD     = TimeStepperTEBD(sites,J,tstep,{"Cutoff=",1E-9});
  OptimalControl<TimeStepperTEBD,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);

  auto times    = generateRange(0,tstep,T);
  auto control  = linspace(cstart,cend,times.size());
  auto Agrad    = OC.getAnalyticGradient(control);
  auto Ngrad    = OC.getNumericGradient(control);

  for (size_t i = 0; i < Agrad.second.size(); i++) {
    std::vector<double> tmp;
    tmp.push_back(Agrad.second.at(i));
    tmp.push_back(Ngrad.second.at(i));
    tmp.push_back( abs(Agrad.second.at(i)-Ngrad.second.at(i)) );

    result.push_back(tmp);
  }

  return result;
}



int main(){
  int N         = 5;
  int Npart     = 5;
  int locDim    = 5;

  double J      = 1.0;
  double cstart = 3.0;
  double cend   = 10.0;
  double T      = 2.0;

  auto sites    = Boson(N,locDim);
  auto psi_i    = InitializeState(sites,Npart,J,cstart);
  auto psi_f    = InitializeState(sites,Npart,J,cend);

  auto tsteps   = linspace(5e-5,5e-4,8);
  auto data     = testCostPlusFidelity(sites,psi_i,psi_f,tsteps,cstart,cend,T,J);
  saveData(data,"tstep_cost_varfidelity2.txt");

  // std::vector<double> tsteps;
  // tsteps.push_back(1e-2);
  //
  //
  // for (size_t order = 0; order <= 2; order++) {
  //   for (auto& dt : tsteps){
  //     auto data = matchGradients(sites,psi_i,psi_f,dt,cstart,cend,T,J,order);
  //     std::string name = "Gradients_order" + std::to_string(order) + "_tstep" + std::to_string(dt) + ".txt";
  //     saveData(data,name);
  //   }
  // }


  return 0;
}
