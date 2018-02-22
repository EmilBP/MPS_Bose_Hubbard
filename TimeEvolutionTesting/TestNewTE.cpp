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

matrix convertToColumn(const matrix& data){
  matrix cdata;
  for (size_t i = 0; i < data.at(0).size(); i++) {
    std::vector<double> tmp;

    for (auto& row : data){
      tmp.push_back(row.at(i));
    }
    cdata.push_back(tmp);
  }
  return cdata;
}

GPdata convertToGPdata(const matrix& data){
  GPdata result;
  auto cdata = convertToColumn(data);
  auto x = cdata.at(0);

  for (size_t i = 1; i < cdata.size(); i++) {
    std::vector< std::pair<double, double> > tmp;
    auto col = cdata.at(i);

    for (size_t j = 0; j < col.size(); j++) {
      tmp.push_back(std::make_pair( x.at(j) , col.at(j) ));
    }
    result.push_back(tmp);
  }
  return result;
}


matrix
testCostPlusFidelity( SiteSet sites,
                      IQMPS psi_i,
                      IQMPS psi_f,
                      std::vector<double> tsteps,
                      double cstart,
                      double cend,
                      double T,
                      double J,
                      bool plot)
{

  matrix result;
  size_t count  = 1;
  auto H_BH     = HamiltonianBH(sites,J,0);

  for (auto& dt: tsteps){
    auto times    = generateRange(0,dt,T);
    auto control  = linspace(cstart,cend,times.size());

    std::cout << "Calculating using control of length " << control.size() << " w. final value " << control.back() << '\n';

    auto TEBD     = TimeStepperTEBDfast(sites,J,dt,{"Cutoff=",1E-8});
    OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);
    auto CF       = OC.checkCostPlusFidelity(control);
    double var    = calculateVariance(CF.second);

    std::vector<double> tmp;
    tmp.push_back(dt);
    tmp.push_back(CF.first);
    tmp.push_back(var);

    result.push_back(tmp);
    std::cout << "Done " << count++ << " out of " << tsteps.size() << '\n';
  }

  if (plot) {
    auto gpdat  = convertToGPdata(result);
    Gnuplot gp3;
    gp3 << "set xlabel 'dt'\n";
    gp3 << "set ylabel 'Cost'\n";
    gp3 << "plot"
    << gp3.file1d(gpdat.at(0)) << "with linespoints ls 1,"
    << std::endl;

    Gnuplot gp4;
    gp4 << "set xlabel 'dt'\n";
    gp4 << "set ylabel 'var(F)'\n";
    gp4 << "plot"
    << gp4.file1d(gpdat.at(1)) << "with linespoints ls 1,"
    << std::endl;
  }

  return result;
}

matrix
matchGradients(       SiteSet sites,
                      IQMPS psi_i,
                      IQMPS psi_f,
                      double tstep,
                      double cstart,
                      double cend,
                      double T,
                      double J,
                      size_t order)
{

  matrix result;
  auto H_BH     = HamiltonianBH(sites,J,tstep,order);
  auto TEBD     = TimeStepperTEBDfast(sites,J,tstep,{"Cutoff=",1E-9});
  OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);

  auto times    = generateRange(0,tstep,T);
  auto control  = linspace(cstart,cend,times.size());
  auto Agrad    = OC.getAnalyticGradient(control);
  // auto Ngrad    = OC.getNumericGradient(control);

  for (size_t i = 0; i < Agrad.second.size(); i++) {
    std::vector<double> tmp;
    tmp.push_back(Agrad.second.at(i));
    // tmp.push_back(Ngrad.second.at(i));
    // tmp.push_back( abs(Agrad.second.at(i)-Ngrad.second.at(i)) );

    result.push_back(tmp);
  }

  return result;
}

matrix
compareTEalgorithms( SiteSet sites,
                      IQMPS psi_i,
                      IQMPS&psi_f,
                      std::vector<double> tsteps,
                      double cstart,
                      double cend,
                      double T,
                      double J)
{

  matrix result;
  auto H_BH     = HamiltonianBH(sites,J,0);

  for (auto& dt : tsteps){

    std::vector<double> tmp;
    tmp.push_back(dt);

    auto times    = generateRange(0,dt,T);
    auto control  = linspace(cstart,cend,times.size());


    auto TEBD1    = TimeStepperTEBDnew(sites,J,dt,{"Cutoff=",1E-8});
    OptimalControl<TimeStepperTEBDnew,HamiltonianBH> OC1(psi_f,psi_i,TEBD1,H_BH, 0);
    std::cout << "Re-ordered Time-Evolution Algorithm\n";
    clock_t begin = clock();
    auto CF       = OC1.checkCostPlusFidelity(control);
    clock_t end   = clock();
    std::cout << "Runtime = " <<  double(end - begin) / CLOCKS_PER_SEC << '\n';
    double var    = calculateVariance(CF.second);
    tmp.push_back(CF.first);
    tmp.push_back(var);
    tmp.push_back((end - begin) / CLOCKS_PER_SEC);


    // auto TEBD2    = TimeStepperTEBD(sites,J,dt,{"Cutoff=",1E-8});
    // OptimalControl<TimeStepperTEBD,HamiltonianBH> OC2(psi_f,psi_i,TEBD2,H_BH, 0);
    // std::cout << "Original Time-Evolution Algorithm\n";
    // begin         = clock();
    // CF            = OC2.checkCostPlusFidelity(control);
    // end           = clock();
    // std::cout << "Runtime = " <<  double(end - begin) / CLOCKS_PER_SEC << '\n';
    // var           = calculateVariance(CF.second);
    // tmp.push_back(CF.first);
    // tmp.push_back(var);
    // tmp.push_back((end - begin) / CLOCKS_PER_SEC);

    auto TEBD3    = TimeStepperTEBDfast(sites,J,dt,{"Cutoff=",1E-8});
    OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC3(psi_f,psi_i,TEBD3,H_BH, 0);
    std::cout << "Fast Time-Evolution Algorithm\n";
    begin         = clock();
    CF            = OC3.checkCostPlusFidelity(control);
    end           = clock();
    std::cout << "Runtime = " <<  double(end - begin) / CLOCKS_PER_SEC << '\n';
    var           = calculateVariance(CF.second);
    tmp.push_back(CF.first);
    tmp.push_back(var);
    tmp.push_back((end - begin) / CLOCKS_PER_SEC);


    result.push_back(tmp);
  }

  return result;
}

matrix
testBackwardsPropagation( SiteSet sites,
                      IQMPS psi_i,
                      IQMPS psi_f,
                      double dt,
                      double cstart,
                      double cend,
                      double T,
                      double J,
                      bool plot)
{

  matrix result;

  auto times    = generateRange(0,dt,T);
  auto control  = linspace(cstart,cend,times.size());
  auto TEBD     = TimeStepperTEBDfast(sites,J,dt,{"Cutoff=",1E-8});
  auto H_BH     = HamiltonianBH(sites,J,0);
  OptimalControl<TimeStepperTEBDfast,HamiltonianBH> OC(psi_f,psi_i,TEBD,H_BH, 0);

  std::vector<IQMPS> forwards;
  std::vector<IQMPS> backwards;

  forwards.push_back(psi_i);
  for (size_t i = 0; i < control.size()-1; i++) {
    TEBD.step(psi_i,control.at(i),control.at(i+1),true);
    forwards.push_back(psi_i);
  }

  backwards.push_back(psi_i);
  for (size_t i = control.size()-1; i > 0; i--) {
    TEBD.step(psi_i,control.at(i),control.at(i-1),false);
    backwards.push_back(psi_i);
  }
  std::reverse(backwards.begin(),backwards.end());


  auto OCres = OC.checkCostPlusFidelity(control);


  for (size_t i = 0; i < forwards.size(); i++) {
    std::vector<double> tmp;
    double re, im;
    overlap(forwards.at(i),backwards.at(i),re,im);
    tmp.push_back(times.at(i));
    tmp.push_back(re);
    tmp.push_back(im);
    tmp.push_back(re*re+im*im);
    tmp.push_back(OCres.second.at(i));

    result.push_back(tmp);
  }

  if (plot) {
    auto gpdat = convertToGPdata(result);
    Gnuplot gp1;
    // gp1 << "set yrange [0.9999:1]\n";
    gp1 << "set xlabel 't'\n";
    gp1 << "set ylabel 'Fidelity'\n";
    gp1 << "plot"
    << gp1.file1d(gpdat.at(2)) << "with lines title 'for-back',"
    << std::endl;

    Gnuplot gp2;
    // gp2 << "set yrange [0.94:0.945]\n";
    gp2 << "set xlabel 't'\n";
    gp2 << "set ylabel 'Fidelity'\n";
    gp2 << "plot"
    << gp2.file1d(gpdat.at(3)) << "with lines title 'psi-chi',"
    << std::endl;
  }

  return result;
}

matrix
compareSpeed(         SiteSet sites,
                      IQMPS psi_i,
                      double dt,
                      double cstart,
                      double cend,
                      double T,
                      double J)
{

  matrix result;

  auto times    = generateRange(0,dt,T);
  auto control  = linspace(cstart,cend,times.size());
  clock_t begin, end;
  IQMPS psi;
  std::vector<double> tmp;
  double TT;

  auto TEBDold  = TimeStepperTEBD(sites,J,dt,{"Cutoff=",1E-8});
  auto TEBDnew  = TimeStepperTEBDnew(sites,J,dt,{"Cutoff=",1E-8});
  auto TEBDfast = TimeStepperTEBDfast(sites,J,dt,{"Cutoff=",1E-8});
  auto tMPO     = TimeStepperMPO(sites,J,dt,{"Cutoff=",1E-8});



  psi = psi_i;
  begin = clock();
  for (size_t i = 0; i < control.size()-1; i++) {
    tMPO.step(psi,control.at(i),control.at(i+1),true);
  }
  end = clock();
  TT = double(end - begin) / CLOCKS_PER_SEC;
  tmp.push_back(TT);
  std::cout << "Runtime MPO = " << TT << '\n';

  psi = psi_i;
  begin = clock();
  for (size_t i = 0; i < control.size()-1; i++) {
    TEBDold.step(psi,control.at(i),control.at(i+1),true);
  }
  end = clock();
  TT = double(end - begin) / CLOCKS_PER_SEC;
  tmp.push_back(TT);
  std::cout << "Runtime TEBD (old) = " << TT << '\n';

  psi = psi_i;
  begin = clock();
  for (size_t i = 0; i < control.size()-1; i++) {
    TEBDnew.step(psi,control.at(i),control.at(i+1),true);
  }
  end = clock();
  TT = double(end - begin) / CLOCKS_PER_SEC;
  tmp.push_back(TT);
  std::cout << "Runtime TEBD (new) = " << TT << '\n';

  psi = psi_i;
  begin = clock();
  for (size_t i = 0; i < control.size()-1; i++) {
    TEBDfast.step(psi,control.at(i),control.at(i+1),true);
  }
  end = clock();
  TT = double(end - begin) / CLOCKS_PER_SEC;
  tmp.push_back(TT);
  std::cout << "Runtime TEBD (fast) = " << TT << '\n';

  result.push_back(tmp);
  return result;
}


matrix
testTimeEvolution( SiteSet sites,
                      int Npart,
                      IQMPS psi_i,
                      double dt,
                      double T,
                      double Uevol,
                      double Jevol,
                      bool plot)
{

  matrix result;

  auto times    = generateRange(0,dt,T);
  std::vector<double> control(times.size(),2.0);
  auto TEBD     = TimeStepperTEBDfast(sites,Jevol,dt,{"Cutoff=",1E-8});

  std::vector<double> condfrac;

  auto lambda1 = correlationTerm(sites,psi_i,"Adag","A");
  condfrac.push_back((double) lambda1/Npart);

  for (size_t i = 0; i < control.size()-1; i++) {
    TEBD.step(psi_i,control.at(i),control.at(i+1),true);

    lambda1 = correlationTerm(sites,psi_i,"Adag","A");
    condfrac.push_back((double) lambda1/Npart);
  }

  for (size_t i = 0; i < condfrac.size(); i++) {
    std::vector<double> tmp;
    tmp.push_back(times.at(i));
    tmp.push_back(condfrac.at(i));

    result.push_back(tmp);
  }

  if (plot) {
    auto gpdat = convertToGPdata(result);
    Gnuplot gp1;
    gp1 << "set xlabel 't'\n";
    gp1 << "set ylabel 'f_c'\n";
    gp1 << "plot"
    << gp1.file1d(gpdat.front()) << "with lines,"
    << std::endl;
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
  auto psi_SF   = SetupSuperfluid(sites,Npart);

  // //
  // //  Test TE algorithm through SF revival
  // //
  // auto dataTE  = testTimeEvolution(sites,Npart,psi_SF,1e-2,10,1,0.03,true);
  //
  // //
  // //  Compare backwards and forwards propagation
  // //
  // auto dataBF  = testBackwardsPropagation(sites,psi_i,psi_f,1e-2,cstart,cend,T,J,true);
  //
  // //
  // //  Compare cost and variance of fidelity for various time steps
  // //
  // auto tsteps   = linspace(1e-3,1e-2,15);
  // auto dataCF   = testCostPlusFidelity(sites,psi_i,psi_f,tsteps,cstart,cend,T,J,true);


  //
  //  Match gradients
  //
  std::vector<double> tsteps;
  tsteps.push_back(1e-2);
  // tsteps.push_back(1e-3);

  for (size_t order = 0; order <= 2; order++) {
    for (auto& dt : tsteps){
      auto data = matchGradients(sites,psi_i,psi_f,dt,cstart,cend,T,J,order);
      std::string name = "Gradients_order" + std::to_string(order) + "_tstep" + std::to_string(dt) + ".txt";
      saveData(data,name);
    }
  }


  // auto tsteps  = linspace(1e-3,1e-2,10);
  // std::vector<double> tsteps = {1e-2};
  // auto data    = compareTEalgorithms(sites,psi_i,psi_f,tsteps,cstart,cend,T,J);
  // printData(data);
  // saveData(data,"compareTEalgorithmsN15.txt");




  // std::vector<int> sizes = {5,6,7,8,9,10};
  // for (auto& s : sizes){
  //
  //   int N         = s;
  //   int Npart     = s;
  //   int locDim    = 5;
  //
  //   auto sites    = Boson(N,locDim);
  //   auto psi_i    = InitializeState(sites,Npart,J,cstart);
  //   auto psi_f    = InitializeState(sites,Npart,J,cend);
  //   auto data     = compareSpeed(sites,psi_i,1e-2,cstart,cend,T,J);
  //   std::string name = "compareTEalgorithmsN" + std::to_string(s) + "d" + std::to_string(locDim) + ".txt";
  //   saveData(data,name);
  // }

  return 0;
}
