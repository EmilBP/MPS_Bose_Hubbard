#include "BoseHubbardMPO.h"

BoseHubbardMPO::BoseHubbardMPO(SiteSet& sites, double overlayTrapDepth)
  : baseMPO(AutoMPO(sites)), N(sites.N()), overlayTrapDepth(overlayTrapDepth), L(sites.N()) {

  hasData = false;
}

void BoseHubbardMPO::interpolateData(double V0, double& U, double& J, double& dUdV, double& dJdV){
  size_t index = 0;

  if (V0 < UJdata[0][0] || V0 > UJdata[UJdata.size()-1][0]) {
    std::cout << "Error: V0 not within data!" << '\n';
    return;
  }

  while(UJdata[index][0] < V0 && index < UJdata.size()) {
    index++;
  }

  dUdV   = (UJdata[index][1]-UJdata[index-1][1])/(UJdata[index][0]-UJdata[index-1][0]);
  dJdV   = (UJdata[index][2]-UJdata[index-1][2])/(UJdata[index][0]-UJdata[index-1][0]);
  U      = UJdata[index-1][1] + dUdV*(V0-UJdata[index-1][0]);
  J      = UJdata[index-1][2] + dJdV*(V0-UJdata[index-1][0]);
}

void BoseHubbardMPO::loadUJdata(std::string& filename){
  UJdata.clear();
  std::ifstream ifs(filename);
  std::string tempstr;
  double var1, var2, var3;

  while (std::getline(ifs, tempstr)) {
    std::istringstream iss(tempstr);
    std::vector<double> tempv;
    while (iss >> var1 >> var2 >> var3) {
      tempv = {var1 , var2, var3};
    }
    UJdata.push_back(tempv); // row based
  }
  hasData = true;
}

AutoMPO BoseHubbardMPO::updateMPO(double control){
  AutoMPO ampo  = baseMPO;
  double U, J, dU, dJ;
  interpolateData(control,U,J,dU,dJ);

  for(int i = 1; i < N; ++i) {
    ampo += -J,"A",i,"Adag",i+1;
    ampo += -J,"Adag",i,"A",i+1;
  }
  for (int i = 1; i <= N; ++i) {
    ampo += U/2.0,"N(N-1)",i;
    ampo += overlayTrapDepth*(i-0.5*(L-1))*(i-0.5*(L-1)),"N",i;
  }

  return ampo;
}

AutoMPO BoseHubbardMPO::dHdu(double control){
  AutoMPO ampo = baseMPO;
  double U, J, dUdV, dJdV;
  interpolateData(control,U,J,dUdV,dJdV);

  for(int i = 1; i < N; ++i) {
    ampo += -dJdV,"A",i,"Adag",i+1;
    ampo += -dJdV,"Adag",i,"A",i+1;
  }
  for (int i = 1; i <= N; ++i) {
    ampo += dUdV/2.0,"N(N-1)",i;
  }

  return ampo;
}
