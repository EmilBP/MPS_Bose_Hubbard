#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

typedef std::vector<std::vector<double>> matrix;

void interpolate(matrix& data, double V0, double& U, double& J){
  size_t index = 0;

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

  double U, J, V0 = 4.5;
  interpolate(v,V0,U,J);
  std::cout << "V0 = " << V0 << ", U = " << U << ", J = " << J << '\n';
}
