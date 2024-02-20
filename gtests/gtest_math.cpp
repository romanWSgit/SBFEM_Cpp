//
// Created by Roman Wallner- Silberhuber on 15.02.24.
//
#include "sbfem_math.h"
#include <Eigen/Dense>
#include <gtest/gtest.h>

class LagrangeTest : public ::testing::Test
{
  protected:
    // You can remove any or all of the following functions if their bodies
    // would be empty.

    LagrangeTest ()
    {
        // You can do set-up work for each test here.
    }

    ~LagrangeTest () override
    {
        // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    void
    SetUp () override
    {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }

    void
    TearDown () override
    {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

    // Class members declared here can be used by all tests in the test suite
    // for Foo.
};

// Testen, ob die Funktion den richtigen Wert zurückgibt,
// wenn positive Zahlen übergeben werden.

TEST_F (LagrangeTest, LagrangeTestSuccessfull)
{
    Eigen::VectorXd xm = evenDistributedLagrangianPoints (3);
    double x = 0.2;
    int i = 0;
    double expected_output = -0.032;
    EXPECT_NEAR (lagrange (x, i, xm), expected_output, 0.001);
}

TEST_F (LagrangeTest, LagrangeTestSuccessfullAt1)
{
    Eigen::VectorXd xm = evenDistributedLagrangianPoints (3);
    double x = 0.2;
    int i = 0;
    double expected_output = -0.032;
    EXPECT_NEAR (lagrange (x, i, xm), expected_output, 0.001);
}
