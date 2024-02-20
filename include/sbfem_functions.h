/**
 * \file   sbfem_functions.h
 * \author Roman Wallner- Silberhuber
 * \date   04.06.2023
 * \brief  This header contains the sbfem transformation of geometry functions
 */

#ifndef SBFEM_SBFEM_FUNCTIONS_H
#define SBFEM_SBFEM_FUNCTIONS_H

#include "ShapeFunctionType.h"
#include <Eigen/Dense>
#include <iostream> // for std::cerr
#include <tuple>    // for std::tie

Eigen::Vector2d r_hat_c(
    double xi, double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
    ShapeFunctionType shape_function,
    const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

Eigen::Vector2d r_hat(double xi, double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function);

Eigen::Vector2d r_c(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                    ShapeFunctionType shape_function_type,
                    const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

Eigen::Vector2d r(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                  ShapeFunctionType shape_function_type,
                  const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

Eigen::Matrix2d j_mat(double eta, const Eigen::VectorXd &coord_vec,
                      int poly_ord, ShapeFunctionType shape_function_type,
                      const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

#endif // SBFEM_SBFEM_FUNCTIONS_H
