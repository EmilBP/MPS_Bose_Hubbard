#include <gtest/gtest.h>
#include "itensor/all.h"

#include "../include/boson.h"
#include "../include/correlations.h"


using namespace itensor;

struct twoSiteCorrelationTest : testing::Test {
  SiteSet* sites;
  IQMPS* psi;

  twoSiteCorrelationTest() {
    int N = 2;
    int d = 5;

    sites = new Boson(N,d);
    auto state = InitState(*sites);
    state.set(1,"Occ1");
    state.set(2,"Occ2");

    psi = new IQMPS(state);
  }

  ~twoSiteCorrelationTest() {
    delete sites;
    delete psi;
  }
};

struct operatorTest : testing::Test {
  SiteSet* sites;
  IQMPS* psi;
  IQMPO* N_opp;
  IQMPO* A_opp;
  IQMPO* Adag_opp;

  operatorTest() {
    int N = 2;
    int d = 5;

    sites = new Boson(N,d);
    auto state = InitState(*sites);
    state.set(1,"Occ2");
    state.set(2,"Occ2");

    psi = new IQMPS(state);

    auto ampoN = AutoMPO(*sites);
    ampoN += 1,"N",1;
    N_opp = new IQMPO(ampoN);

    auto ampoAdag = AutoMPO(*sites);
    ampoAdag += 1,"Adag",1,"A",2;
    Adag_opp = new IQMPO(ampoAdag);

    auto ampoA = AutoMPO(*sites);
    ampoA += 1,"A",1,"Adag",2;
    A_opp = new IQMPO(ampoA);

  }

  ~operatorTest() {
    delete sites;
    delete psi;
    delete N_opp;
    delete A_opp;
    delete Adag_opp;
  }
};

struct DISABLED_GroundStateEnergyTest : testing::Test {
  SiteSet* sites;
  IQMPS* psi;
  int N;
  int d;

  DISABLED_GroundStateEnergyTest() {}

  void Updateparams(int number) {
    N=number;
    d = N;
    sites = new Boson(N,d);
    auto state = InitState(*sites);
    //state.set(1,"Occ1");
    //state.set(2,"Occ2");

    int p = N;
    for(int i = N; i >= 1; --i)
        {
        if (p >= 1) {
          //println("Singly occupying site ",i);
          state.set(i,"Occ1");
          --p;
        }
        else
            {
            state.set(i,"Emp");
            }
        }

    psi = new IQMPS(state);
  }

  ~DISABLED_GroundStateEnergyTest() {
    delete sites;
    delete psi;
  }

  double updateHamiltonian(double J, double U, double eps=0) {
    auto ampo = AutoMPO(*sites);
    for(int i = 1; i < N; ++i) {
      ampo += J,"A",i,"Adag",i+1;
      ampo += J,"Adag",i,"A",i+1;
    }
    for (int i = 1; i <= N; ++i) {
      ampo += U/2.0,"N(N-1)",i;
      ampo += eps,"N",i;
    }

    auto H = IQMPO(ampo);

    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,50,50,100,200;
    sweeps.cutoff() = 1E-9;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;

    auto energy = dmrg(*psi,H,sweeps,"Quiet");
    return energy;
  }
};

struct singleParticleDensityMatrixTest : testing::Test {
  SiteSet* sites;
  IQMPS* psi;

  singleParticleDensityMatrixTest(){ }
  ~singleParticleDensityMatrixTest() {
    delete sites;
    delete psi;
  }

  void setupMott(int fillingFraction, int N){
    int Npart = fillingFraction*N;
    int d = fillingFraction*2;

    sites = new Boson(N,d);
    auto state = InitState(*sites);

    for (size_t i = 1; i <= N; i++) {
        state.set(i,nameint("Occ",fillingFraction));
    }

    psi = new IQMPS(state);

  }

};

TEST_F(operatorTest, testParticleNumber){
  auto args = Args("Cutoff=",1E-9,"Maxm=",10);

  ASSERT_DOUBLE_EQ( 2.0, overlap(*psi,*N_opp,*psi));

  // raise particle number by 1, EVAL N, lower by 1
  *psi = exactApplyMPO(*Adag_opp,*psi,args);
  normalize(*psi);
  ASSERT_DOUBLE_EQ( 3.0, overlap(*psi,*N_opp,*psi));
  *psi = exactApplyMPO(*A_opp,*psi,args);
  normalize(*psi);

  // lower particle number by 1, EVAL N, raise by 1
  *psi = exactApplyMPO(*A_opp,*psi,args);
  normalize(*psi);
  ASSERT_DOUBLE_EQ( 1.0, overlap(*psi,*N_opp,*psi));
  *psi = exactApplyMPO(*Adag_opp,*psi,args);
  normalize(*psi);
}


TEST_F(DISABLED_GroundStateEnergyTest, test5ParticleEnergy) {
    Updateparams(5);
    EXPECT_NEAR(-8.66025, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.319776, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-0.6382, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(DISABLED_GroundStateEnergyTest, test4ParticleEnergy) {
    Updateparams(4);
    EXPECT_NEAR(-6.47214, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.239808, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-0.47846, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(DISABLED_GroundStateEnergyTest, test6ParticleEnergy) {
    Updateparams(6);
    EXPECT_NEAR(-10.8116, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.399744, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-0.797941, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(DISABLED_GroundStateEnergyTest, test7ParticleEnergy) {
    Updateparams(7);
    EXPECT_NEAR(-12.9343, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.479712, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-0.957682, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(DISABLED_GroundStateEnergyTest, test8ParticleEnergy) {
    Updateparams(8);
    EXPECT_NEAR(-15.0351, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.559679, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-1.11742, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(DISABLED_GroundStateEnergyTest, test9ParticleEnergy) {
    Updateparams(9);
    EXPECT_NEAR(-17.119, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.639647, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-1.27716, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(DISABLED_GroundStateEnergyTest, test10ParticleEnergy) {
    Updateparams(10);
    EXPECT_NEAR(-19.1899, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.719615, updateHamiltonian(-1.0, 50), 1.5e-5);
    EXPECT_NEAR(-1.43691, updateHamiltonian(-1.0, 25), 1.5e-5);
}

TEST_F(twoSiteCorrelationTest, testBosonSiteWithOcc1) {
    int site1 = 1;
    int site2 = 1;

    ASSERT_EQ( (Cplx) 1, correlations::correlationFunction(*sites,*psi,"N",site1,"N",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"I",site1,"N(N-1)",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"N(N-1)",site1,"I",site2));

    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"Adag",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"A",site1,"A",site2));
    ASSERT_EQ( (Cplx) 1, correlations::correlationFunction(*sites,*psi,"Adag",site1,"A",site2));
    ASSERT_EQ( (Cplx) 2, correlations::correlationFunction(*sites,*psi,"A",site1,"Adag",site2));

    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"N",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"A",site1,"N",site2));
}

TEST_F(twoSiteCorrelationTest, testBosonSiteWithOcc2) {
    int site1 = 2;
    int site2 = 2;

    ASSERT_EQ( (Cplx) 4.0, correlations::correlationFunction(*sites,*psi,"N",site1,"N",site2));
    ASSERT_EQ( (Cplx) 2.0, correlations::correlationFunction(*sites,*psi,"I",site1,"N(N-1)",site2));
    ASSERT_EQ( (Cplx) 2.0, correlations::correlationFunction(*sites,*psi,"N(N-1)",site1,"I",site2));

    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"Adag",site2));
    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"A",site2));
    ASSERT_EQ( (Cplx) 2.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"A",site2));
    ASSERT_DOUBLE_EQ( 3, real(correlations::correlationFunction(*sites,*psi,"A",site1,"Adag",site2)));

    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"N",site2));
    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"N",site2));
}

TEST_F(twoSiteCorrelationTest, testBosonSitesWithOcc1andOcc2) {
    int site1 = 1;
    int site2 = 2;

    ASSERT_EQ( (Cplx) 2, correlations::correlationFunction(*sites,*psi,"N",site1,"N",site2));
    ASSERT_EQ( (Cplx) 2, correlations::correlationFunction(*sites,*psi,"I",site1,"N(N-1)",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"N(N-1)",site1,"I",site2));

    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"Adag",site2));
    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"A",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"A",site2));
    ASSERT_EQ( (Cplx) 0, correlations::correlationFunction(*sites,*psi,"A",site1,"Adag",site2));

    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"N",site2));
    ASSERT_EQ( (Cplx) 0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"N",site2));
}

TEST_F(singleParticleDensityMatrixTest, testMott){
  for (size_t i = 1; i < 5; i++) {
    for (size_t j = 5; j < 25; j=j+5) {
      setupMott(i,j);
      ASSERT_DOUBLE_EQ(i, correlations::correlationTerm(*sites,*psi,"Adag","A"));
    }
  }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
