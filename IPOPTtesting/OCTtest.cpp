#include "OCdummy_nlp.hpp"
#include "OptimalControlDummy.hpp"
#include "IpIpoptApplication.hpp"

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

void printData(const matrix& data){
  for (auto& row : data){
    for (auto& val : row){
      std::cout << val << "\t";
    }
    std::cout << "\n";
  }
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


int main(){

  double tstep  = 1e-2;
  double T      = 5;
  double cstart = 2;
  double cend   = 7;

  std::vector<double> weights = {5.5 , -1.2 , -6.3 , 0.3};
  auto OCD      = OptimalControlDummy(weights,tstep);

  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<TNLP> mynlp = new OCdummy_nlp(OCD,cstart,cend);

  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-9);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
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
