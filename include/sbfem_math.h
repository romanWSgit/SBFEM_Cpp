//
// Created by Roman Wallner- Silberhuber on 20.05.23.
//

#ifndef SBFEM_SBFEM_MATH_H
#define SBFEM_SBFEM_MATH_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath> // Für std::pow
#include <iostream>
#include <tuple>
#include <string>
#include "ShapeFunctionType.h"

/**
 * @brief Calculate the Lagrange polynomial at a given point.
 *
 * This function computes the value of the Lagrange polynomial for a given
 * point `x` and a specific index `i` within the set of sample points `xm`.
 * The Lagrange polynomial is used for polynomial interpolation where the
 * polynomial of degree `n` passes through `(n+1)` given points.
 *
 * @param x The point at which the Lagrange polynomial is evaluated.
 * @param i The index of the current sample point in `xm` used for calculating
 *          the polynomial. It must be less than the size of `xm`.
 * @param xm A vector (Eigen::VectorXd) containing the sample points. These
 *           points are the `x` coordinates of known data points.
 * @return double The value of the Lagrange polynomial at point `x` for the
 *                given index `i`.
 * @note The size of `xm` should be greater than 1, and `i` should be a valid
 *       index within the bounds of `xm`.
 */
double
lagrange (double x, int i, Eigen::VectorXd xm);

/**
 * @brief Computes even distributed Lagrangian points for a given polynomial
 * order.
 *
 * @param polyOrd The order of the polynomial.
 * @return Eigen::VectorXd A vector of points distributed evenly from -1 to 1,
 * inclusive.
 */
Eigen::VectorXd
evenDistributedLagrangianPoints (int polyOrd);

/**
 * Computes the differentiated Lagrange polynomial evaluated at a given point.
 *
 * @param x The point at which to evaluate the polynomial.
 * @param j The index of the x-coordinate to exclude from the interpolation.
 * @param xm The vector of x-coordinates used for interpolation.
 * @return The value of the differentiated Lagrange polynomial at the given
 * point.
 */
double
lagrange_diff (double x, int j, Eigen::VectorXd xm);

/**
 * @brief Calculate the integrated value of the Legendre polynomial of degree p
 * at x.
 *
 * The Legendre polynomial is calculated based on the degree p and the value of
 * x.
 *
 * @param p The degree of the Legendre polynomial.
 * @param x The value of x at which to evaluate the Legendre polynomial.
 * @return The integrated value of the Legendre polynomial at x.
 */
double
legendre_poly_integrated (double x, int p);

/**
 * Calculates the integrated and then differentiated Legendre polynomials.
 *
 * @param p The order of the Legendre polynomial.
 * @param x The input value.
 * @return The integrated and then differentiated Legendre polynomial at the
 * given order and input value.
 */
double
legendre_poly_integrated_diff (double x, int p);

/**
 * @brief Calculates the shape functions and their gradients at a given
 * point on an element.
 *
 * This function calculates the shape functions and their gradients at a
 * given evaluation point on an element. The shape functions are calculated
 * using Lagrange interpolation with even distributed Lagrangian points. The
 * order of the polynomial basis functions used in the calculation can be
 * specified.
 *
 * @param eta       The evaluation point on the element
 * @param poly_ord  The order of the polynomial basis functions
 * @return          A tuple containing the shape function vector and the
 * gradient matrix
 */
std::tuple<Eigen::VectorXd, Eigen::MatrixXd> shape_N(
    double eta, int poly_ord, ShapeFunctionType shapeFct);

/**
 * @brief Calculates the shape function and its derivative matrix for a given
 * eta value and polynomial order.
 *
 * This function calculates the shape function and its derivative matrix for a
 * given eta value and polynomial order.
 *
 * @param eta The input eta value.
 * @param poly_ord The polynomial order.
 * @return std::tuple<Eigen::VectorXd, Eigen::MatrixXd> A tuple containing the
 * shape function vector and the derivative matrix.
 *         - The shape function vector is a column vector containing the shape
 * function values at the given eta value.
 *         - The derivative matrix is a 2xN matrix (N is the polynomial order)
 * where each column represents the derivative values at the given eta value.
 */
std::tuple<Eigen::VectorXd, Eigen::MatrixXd> shape_dN(
    double eta, int poly_ord, ShapeFunctionType shapeFct);

#endif // SBFEM_SBFEM_MATH_H
