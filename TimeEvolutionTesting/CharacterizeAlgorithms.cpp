#include "itensor/all.h"
#include "boson.h"
#include "HamiltonianBH.hpp"
#include "TimeStepperTEBD.hpp"
#include "TimeStepperTEBDfast.hpp"
#include "TimeStepperTEBDnew.hpp"
#include "TimeStepperMPO.hpp"
#include "InitializeState.hpp"
#include "correlations.hpp"
#include "OptimalControl.hpp"
#include "OptimalControlDummy.hpp"
#include <fstream>
#include <string>
#include <time.h>
#include "gnuplot-iostream.h"

using namespace itensor;

using matrix = std::vector< std::vector<double> >;
using GPdata = std::vector< std::vector< std::pair< double,double> > >;

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

void printData(const matrix& data){
  for (auto& row : data){
    for (auto& val : row){
      std::cout << val << "\t";
    }
    std::cout << "\n";
  }
}

template <class T>
double avgRuntime(T timestepper, int Nruns, std::vector<double> control, IQMPS psi)
{
  std::vector<double> runtimes;
  for (size_t j = 0; j < Nruns; j++) {
    auto psitmp = psi;
    clock_t begin = clock();
    for (size_t i = 0; i < control.size()-1; i++) {
      timestepper.step(psitmp,control.at(i),control.at(i+1),true);
    }
    clock_t end = clock();
    runtimes.push_back(double(end - begin) / CLOCKS_PER_SEC);
    std::cout << "Finished run " << j+1 << '\n';
  }

  double sum = 0;
  for (auto& val : runtimes){
    sum += val;
  }

  return (double) sum/Nruns;
}

std::vector<double> compareSpeed(SiteSet sites, double dt, double T, int Nruns)
{
  std::vector<double> results;

  double cstart = 2.0;
  double cend   = 10.0;
  double J      = 1.0;
  int Npart     = sites.N();

  auto times    = generateRange(0,dt,T);
  auto control  = linspace(cstart,cend,times.size());
  auto psi_i    = InitializeState(sites,Npart,J,cstart);

  auto TEBDold  = TimeStepperTEBD(sites,J,dt,{"Cutoff=",1E-8});
  // auto TEBDnew  = TimeStepperTEBDnew(sites,J,dt,{"Cutoff=",1E-8});
  auto TEBDfast = TimeStepperTEBDfast(sites,J,dt,{"Cutoff=",1E-8});
  auto tMPO     = TimeStepperMPO(sites,J,dt,{"Cutoff=",1E-8});

  results.emplace_back( avgRuntime(TEBDold, Nruns, control, psi_i) );
  // results.emplace_back( avgRuntime(TEBDnew, Nruns, control, psi_i) );
  results.emplace_back( avgRuntime(TEBDfast, Nruns, control, psi_i) );
  results.emplace_back( avgRuntime(tMPO, Nruns, control, psi_i) );

  return results;
}


int main(){

  std::vector<int> Nvec = {5,6,7,8,9,10,11};
  double T    = 1.0;
  double dt   = 1e-2;
  int Nruns   = 5;

  matrix runtimeData1;
  matrix runtimeData2;

  for (auto& N : Nvec){
    int Npart   = N;
    int d1      = N;
    int d2      = Nvec.front();

    auto sites1 = Boson(N,d1);
    auto sites2 = Boson(N,d2);

    runtimeData1.emplace_back( compareSpeed(sites1,dt,T,Nruns) );
    runtimeData2.emplace_back( compareSpeed(sites2,dt,T,Nruns) );

    std::cout << "Finished all runs for N = " << N << '\n';
  }

  saveData(runtimeData1,"RuntimeDataFullDim.txt");
  saveData(runtimeData2,"RuntimeDataFixedDim.txt");

  return 0;
}
