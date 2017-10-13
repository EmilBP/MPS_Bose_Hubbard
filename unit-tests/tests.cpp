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
