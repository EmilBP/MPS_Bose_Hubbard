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

    sites = new Boson(N);
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

struct GroundStateEnergyTest : testing::Test {
  SiteSet* sites;
  IQMPS* psi;
  int N;

  GroundStateEnergyTest() {}

  void Updateparams(int number) {
    N=number;
    sites = new Boson(N);
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

  ~GroundStateEnergyTest() {
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
      ampo += U/2.0,"N",i,"N-1",i;
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

TEST_F(GroundStateEnergyTest, test5ParticleEnergy) {
    Updateparams(5);
    EXPECT_NEAR(-8.66025, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.319776, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-0.6382, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(GroundStateEnergyTest, test4ParticleEnergy) {
    Updateparams(4);
    EXPECT_NEAR(-6.47214, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.239808, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-0.47846, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(GroundStateEnergyTest, test6ParticleEnergy) {
    Updateparams(6);
    EXPECT_NEAR(-10.8116, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.399744, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-0.797941, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(GroundStateEnergyTest, test7ParticleEnergy) {
    Updateparams(7);
    EXPECT_NEAR(-12.9343, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.479712, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-0.957682, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(GroundStateEnergyTest, test8ParticleEnergy) {
    Updateparams(8);
    EXPECT_NEAR(-15.0351, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.559679, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-1.11742, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(GroundStateEnergyTest, test9ParticleEnergy) {
    Updateparams(9);
    EXPECT_NEAR(-17.119, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.639647, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-1.27716, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(GroundStateEnergyTest, test10ParticleEnergy) {
    Updateparams(10);
    EXPECT_NEAR(-19.1899, updateHamiltonian(-1.0, 0), 1.5e-5);
    EXPECT_NEAR(-0.719615, updateHamiltonian(-1.0, 50), 1.5e-5);    
    EXPECT_NEAR(-1.43691, updateHamiltonian(-1.0, 25), 1.5e-5);    
}

TEST_F(twoSiteCorrelationTest, testBosonSiteWithOcc1) {
    int site1 = 1;
    int site2 = 1;

    ASSERT_DOUBLE_EQ(1, correlations::correlationFunction(*sites,*psi,"N",site1,"N",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"N",site1,"N-1",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"N-1",site1,"N",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"N-1",site1,"N-1",site2));

    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"Adag",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"A",site1,"A",site2));
    ASSERT_DOUBLE_EQ(1, correlations::correlationFunction(*sites,*psi,"Adag",site1,"A",site2));
    ASSERT_DOUBLE_EQ(2, correlations::correlationFunction(*sites,*psi,"A",site1,"Adag",site2));

    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"N",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"A",site1,"N",site2));
}

TEST_F(twoSiteCorrelationTest, testBosonSiteWithOcc2) {
    int site1 = 2;
    int site2 = 2;

    ASSERT_DOUBLE_EQ(4.0, correlations::correlationFunction(*sites,*psi,"N",site1,"N",site2));
    ASSERT_DOUBLE_EQ(2.0, correlations::correlationFunction(*sites,*psi,"N",site1,"N-1",site2));
    ASSERT_DOUBLE_EQ(2.0, correlations::correlationFunction(*sites,*psi,"N-1",site1,"N",site2));
    ASSERT_DOUBLE_EQ(1.0, correlations::correlationFunction(*sites,*psi,"N-1",site1,"N-1",site2));

    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"Adag",site2));
    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"A",site2));
    ASSERT_DOUBLE_EQ(2.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"A",site2));
    ASSERT_DOUBLE_EQ(3.0, correlations::correlationFunction(*sites,*psi,"A",site1,"Adag",site2));

    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"N",site2));
    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"N",site2));
}

TEST_F(twoSiteCorrelationTest, testBosonSitesWithOcc1andOcc2) {
    int site1 = 1;
    int site2 = 2;

    ASSERT_DOUBLE_EQ(2, correlations::correlationFunction(*sites,*psi,"N",site1,"N",site2));
    ASSERT_DOUBLE_EQ(1, correlations::correlationFunction(*sites,*psi,"N",site1,"N-1",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"N-1",site1,"N",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"N-1",site1,"N-1",site2));

    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"Adag",site2));
    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"A",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"A",site2));
    ASSERT_DOUBLE_EQ(0, correlations::correlationFunction(*sites,*psi,"A",site1,"Adag",site2));

    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"Adag",site1,"N",site2));
    ASSERT_DOUBLE_EQ(0.0, correlations::correlationFunction(*sites,*psi,"A",site1,"N",site2));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
