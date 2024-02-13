//
// Created by Roman Wallner- Silberhuber on 26.11.23.
//

#include "sbfem_math.h"
#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE (MathSuite1)

Eigen::VectorXd
GetSampleInterpolationPointsVector ()
{
    Eigen::VectorXd xm = evenDistributedLagrangianPoints (3);
    return xm;
}

void
TestLagrangeFunction (Eigen::VectorXd &xm, double x, int i,
                      double expected_output)
{
    // Call the function and check the result
    double result = lagrange (x, i, xm);
    BOOST_CHECK_CLOSE (result, expected_output, 0.001);
}

void
TestShapeFunctions (Eigen::VectorXd &xm, double x, int i,
                    double expected_output)
{
    // Call the function and check the result
    double result = lagrange (x, i, xm);
    BOOST_CHECK_CLOSE (result, expected_output, 0.001);
}

BOOST_AUTO_TEST_CASE (test_lagrange_function)
{
    Eigen::VectorXd xm = GetSampleInterpolationPointsVector ();
    double x = 0.2;
    int i = 0;
    double expected_output = -0.032;
    TestLagrangeFunction (xm, x, i, expected_output);
}

BOOST_AUTO_TEST_CASE (test__shape_functions)
{
    Eigen::VectorXd xm = GetSampleInterpolationPointsVector ();
    double x = 0.2;
    int i = 0;
    double expected_output = -0.032;
    TestLagrangeFunction (xm, x, i, expected_output);
}
BOOST_AUTO_TEST_SUITE_END ()