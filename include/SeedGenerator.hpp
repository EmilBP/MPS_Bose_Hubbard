#ifndef SEEDGENERATOR_HPP
#define SEEDGENERATOR_HPP

#include <vector>
#include <assert.h>
#include "math.h"
#include <stdlib.h>
#include <time.h>


class SeedGenerator{

private:
  static double randomDouble(double min, double max);

public:

  static std::vector<double> linspace(double a, double b, int n);
  static std::vector<double> generateRange(double a, double b, double c);
  static std::vector<double> sigmoid(std::vector<double>& x, double k, double offset);
  static std::vector<double> linsigmoidSeed(double u_start, double u_end, size_t length);
  static std::vector<double> adiabaticSeed(double u_start, double u_end, size_t length);
};

std::vector<double> SeedGenerator::linspace(double a, double b, int n){
  std::vector<double> array;
  double step = (b-a) / (n-1);

  while(a <= b + 1e-7) {
      array.push_back(a);
      a += step;           // could recode to better handle rounding errors
  }
  return array;
}

std::vector<double> SeedGenerator::generateRange(double a, double b, double c) { //equiv to a:b:c
    std::vector<double> array;
    while(a <= c + 1e-7) {
        array.push_back(a);
        a += b;         // could recode to better handle rounding errors
    }
    return array;
}

std::vector<double> SeedGenerator::sigmoid(std::vector<double>& x, double k, double offset){
  std::vector<double> S;
  for (auto& xval : x){
    S.push_back(1.0/(1+exp(-k*(xval-offset))));
  }
  return S;
}

double SeedGenerator::randomDouble(double min, double max){
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}

std::vector<double> SeedGenerator::linsigmoidSeed(double u_start, double u_end, size_t length){
  // f = a*x + b/(1+exp(-c*(x-d))) + e
  auto x  = linspace(0,100,length);
  auto a  = randomDouble(0.01,0.15);
  auto b  = u_end-u_start-a*x.back();
  auto c  = randomDouble(0.06,0.18);
  auto d  = randomDouble(60,80);

  auto S1 = sigmoid(x,0.7,5);
  auto S2 = sigmoid(x,-0.9,100-7);

  for (size_t i = S1.size()/2; i < S1.size(); i++) {
    S1.at(i) = S2.at(i);
  }
  S1.front() = 0;
  S1.back()  = 0;

  size_t i = 0;
  for (auto& fx : x){
    fx = S1.at(i)*(a*fx + b/(1+exp(-c*(fx-d))) + u_start) + (1-S1.at(i))*( (u_end-u_start)/(1+exp(-0.2*(fx-40))) + u_start);
    i++;
  }

  return x;
}

std::vector<double> SeedGenerator::adiabaticSeed(double u_start, double u_end, size_t length){
  auto xlist  = linspace(0,100,length);
  double p    = 3.5;
  double k    = 1.0/3.0;
  double xs   = 40;
  double a    = 0.01;

  for (auto& x : xlist){
    if ( x < xs){
      x = (p - u_start -a*xs)/(1+exp( -k * (x - xs/2.0 ))) + u_start + a*x;
    } else {
      x = exp(log(u_end - p + 1)/(100-xs) *(x-xs)) + p -1;
    }
  }

  return xlist;
}


#endif
