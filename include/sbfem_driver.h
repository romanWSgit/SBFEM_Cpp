//
// Created by Roman Wallner- Silberhuber on 20.02.24.
//

#ifndef SBFEM_DRIVER_H
#define SBFEM_DRIVER_H

#include "MaterialData.h"
#include "SuperElementJson.h"
#include "helper_functions.h"
#include "sbfem_functions.h"
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <functional>
#include <iostream>
#include <omp.h>

Eigen::VectorXd exampleFunction(double x);

double f_new(double t);

template <typename Func>
double integrateScalarKronrod(Func f, double a, double b)
{
    double error_estimate;
    auto q = boost::math::quadrature::gauss_kronrod<double, 15>();
    double result = q.integrate(f, a, b, 5, 1e-10, &error_estimate);
    std::cout << " integration error estimate: " << error_estimate << std::endl;
    return result;
}
template <typename Func2>
double integrateScalarKronrod2D(Func2 f, double ax, double bx, double ay,
                                double by)
{
    auto fx = [&](double x) {
        // Integrate over y for a given x
        auto fy = [&](double y) { return f(x, y); };
        return integrateScalarKronrod(fy, ay, by);
    };
    // Integrate over x
    return integrateScalarKronrod(fx, ax, bx);
}

double integrateComponent(double a, double b, int componentIndex);
double integrateComponentKronrod(double a, double b, int componentIndex);

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, double, double> matrix_integrator(
    double a, double b, const std::function<Eigen::MatrixXd(double)> &matFunc,
    int rows, int cols);

/**
 * @brief Main loop function that performs calculations for each super element.
 *
 * @param sE The SuperElementJson object containing the super element data.
 * @param material_data_list The list of material data to be used for
 * calculations.
 */
void main_loop(const SuperElementJson &sE,
               std::vector<StructData> &material_data_list);

#endif // SBFEM_DRIVER_H