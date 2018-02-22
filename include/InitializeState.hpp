#ifndef INITIALIZESTATE_H
#define INITIALIZESTATE_H

#include "itensor/all.h"
#include "boson.h"

namespace itensor{

inline IQMPS SetupSuperfluid(SiteSet const& sites, int Npart){
  int N = sites.N();

  if (Npart > N) { printf("Npart > N not supported\n"); }

  auto state = InitState(sites);
  int p = Npart;
  for(int i = N; i >= 1; --i)
      {
      if (p >= 1) {
        state.set(i,"Occ1");
        p -= 1;
      }
      else
          {
          state.set(i,"Emp");
          }
      }
  auto psi = IQMPS(state);

  double J = -1.0;
  auto ampo = AutoMPO(sites);
  for(int i = 1; i < N; ++i) {
    ampo += J,"A",i,"Adag",i+1;
    ampo += J,"Adag",i,"A",i+1;
  }
  auto H = IQMPO(ampo);

  auto sweeps = Sweeps(10);
  sweeps.maxm() = 10,20,50,100,200,250,300;
  sweeps.cutoff() = 1E-11;
  sweeps.niter() = 2;
  sweeps.noise() = 1E-7,1E-8,0.0;
  println(sweeps);

  auto energy = dmrg(psi,H,sweeps,{"Quiet",true});
  printf("Energy w.r.t. superfluid Hamiltonian: %f\n",energy);

  return psi;
}

inline IQMPS SetupMottInsulator(SiteSet const& sites, int Npart){
  int N = sites.N();

  if (Npart > N) { printf("Npart > N not supported\n"); }

  auto state = InitState(sites);
  int p = Npart;
  for(int i = N; i >= 1; --i)
      {
      if (p >= 1) {
        state.set(i,"Occ1");
        p -= 1;
      }
      else
          {
          state.set(i,"Emp");
          }
      }
  auto psi = IQMPS(state);

  auto ampo = AutoMPO(sites);
  for(int i = 1; i < N; ++i) {
    ampo += -1e-3,"A",i,"Adag",i+1;
    ampo += -1e-3,"Adag",i,"A",i+1;
  }
  for(int i = 1; i <= N; ++i) {
    ampo += 0.5,"N(N-1)",i;
  }
  auto H = IQMPO(ampo);

  auto sweeps = Sweeps(10);
  sweeps.maxm() = 10,20,50,100,200;
  sweeps.cutoff() = 1E-9;
  sweeps.niter() = 2;
  sweeps.noise() = 1E-7,1E-8,0.0;
  println(sweeps);

  auto energy = dmrg(psi,H,sweeps,{"Quiet",true});
  printf("Energy w.r.t. Mott-Insulator Hamiltonian: %f\n",energy);

  return psi;
}

inline IQMPS InitializeState(const SiteSet& sites, const int Npart, const double J, const double U){
  int N = sites.N();
  auto state = InitState(sites);
  int p = Npart;

  if (Npart > N) { printf("Npart > N not supported\n"); }
  for(int i = N; i >= 1; --i)
      {
      if (p >= 1) {
        state.set(i,"Occ1");
        p -= 1;
      }
      else
          {
          state.set(i,"Emp");
          }
      }
  auto psi = IQMPS(state);

  auto ampo = AutoMPO(sites);
  for(int i = 1; i < N; ++i) {
    ampo += -J,"A",i,"Adag",i+1;
    ampo += -J,"Adag",i,"A",i+1;
  }
  for(int i = 1; i <= N; ++i) {
    ampo += 0.5*U,"N(N-1)",i;
  }
  auto H = IQMPO(ampo);

  auto sweeps = Sweeps(10);
  sweeps.maxm() = 10,20,50,100,200;
  sweeps.cutoff() = 1E-9;
  sweeps.niter() = 2;
  sweeps.noise() = 1E-7,1E-8,0.0;
  
  auto energy = dmrg(psi,H,sweeps,{"Quiet",true});
  return psi;
}

} // end namespace
#endif
