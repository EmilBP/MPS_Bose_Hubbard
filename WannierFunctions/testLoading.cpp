#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>

typedef std::vector<std::vector<double>> matrix;

void interpolate(matrix& data, double V0, double& U, double& J){
  size_t index = 0;

  if (V0 < data[0][0] || V0 > data[data.size()-1][0]) {
    std::cout << "Error: V0 not within data!" << '\n';
    return;
  }

  while(data[index][0] < V0 && index < data.size()) {
    index++;
  }

  double dUdV = (data[index][1]-data[index-1][1])/(data[index][0]-data[index-1][0]);
  double dJdV = (data[index][2]-data[index-1][2])/(data[index][0]-data[index-1][0]);
  U           = data[index-1][1] + dUdV*(V0-data[index-1][0]);
  J           = data[index-1][2] + dJdV*(V0-data[index-1][0]);
}

int main() {
  matrix v;
  std::ifstream ifs("UJparams");
  std::string tempstr;
  double var1, var2, var3;

  while (std::getline(ifs, tempstr)) {
    std::istringstream iss(tempstr);
    std::vector<double> tempv;
    while (iss >> var1 >> var2 >> var3) {
      tempv = {var1 , var2, var3};
    }
    v.push_back(tempv);
  }

  std::vector<double> V0list;
  std::vector<double> Ulist;
  std::vector<double> Jlist;
  for (size_t i = 0; i < v.size(); i++) {
    V0list.push_back(v[i][0]);
    Ulist.push_back(v[i][1]);
    Jlist.push_back(v[i][2]);
  }

  Gnuplot gp;

  // Data will be sent via a temporary file.  These are erased when you call
  // gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
  // (e.g. "gp.file1d(pts, 'mydata.dat')"), then the named file will be created
  // and won't be deleted (this is useful when creating a script).
  gp << "plot '-' with lines title 'U', '-' with lines title 'J'\n";
  gp.send1d(boost::make_tuple(V0list,Ulist));
  gp.send1d(boost::make_tuple(V0list,Jlist));
}
