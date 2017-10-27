#include "correlations.h"
#include <iostream>

Cplx correlations::correlationFunction(SiteSet const& sites, IQMPS& psi, std::string const& opname1, int i, std::string const& opname2, int j){
  auto op_i = sites.op(opname1,i);
  auto op_j = sites.op(opname2,j);


  if (j == i) {
    psi.position(i);

    auto ket = psi.A(i)* op_j*prime(op_i,Site);
    auto bra = prime(prime(dag(psi.A(i)),Site),Site);
    return (bra*ket).real();
  }
  else
  if (j > i) {   }
  else
  if (i > j) {
    int tmp = i;
    i = j;
    j = tmp;
  }

  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  psi.position(i);

  //psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

  //index linking i to i+1:
  auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);

  auto C = psi.A(i)*op_i*prime(dag(psi.A(i)),Site,ir);
  for(int k = i+1; k < j; ++k){
    C *= psi.A(k);
    C *= prime(dag(psi.A(k)),Link);
  }
  C *= psi.A(j);
  C *= op_j;
  //index linking j to j-1:
  auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
  C *= prime(dag(psi.A(j)),jl,Site);

  return C.cplx(); //or C.cplx() if expecting complex
}

ITensor correlations::correlationMatrix(SiteSet const& sites, IQMPS& psi, std::string const& opname1, std::string const& opname2){
  int N = sites.N();
  Index rho_i("rho i",N),
        rho_j = prime(rho_i);

  ITensor rho(rho_i,rho_j);
  Cplx Cij;
  for (int i = 1; i <= N; ++i) {
    Cij = correlationFunction(sites,psi,opname1,i,opname2,i);
    rho.set(rho_i(i),rho_j(i), Cij);
    for (int j = i+1; j <= N; ++j) {
      Cij = correlationFunction(sites,psi,opname1,i,opname2,j);
      rho.set(rho_i(i),rho_j(j), Cij);
      rho.set(rho_i(j),rho_j(i), conj(Cij));
    }
  }
  return rho;
}

Real correlations::correlationTerm(SiteSet const& sites, IQMPS & psi, std::string const& opname1, std::string const& opname2){
  auto rho = correlationMatrix(sites,psi,opname1,opname2);

  ITensor V, D;
  diagHermitian(rho,V,D, {"Maxm",1});

  auto indices = D.inds();
  auto index1 = indices.index(1);
  auto index2 = indices.index(2);

  return D.real(index1(1),index2(1));
}
