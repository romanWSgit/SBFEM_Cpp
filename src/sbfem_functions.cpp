/**
 * \file   sbfem_functions.cpp
 * \author Roman Wallner- Silberhuber
 * \date   04.06.2023
 * \brief  This file contains the Sbfem - transformation of geometry functions
 */

#include "sbfem_functions.h"
#include "sbfem_math.h"

Eigen::Vector2d
r_hat_c (double xi, double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
         const Eigen::Vector2d &centre, std::string_view shape_function)
{
    Eigen::VectorXd shape_N_res;
    Eigen::MatrixXd unused;
    std::tie (shape_N_res, unused) = shape_N (shape_function, eta, poly_ord);
    return xi * shape_N_res[1] * coord_vec + centre;
}

Eigen::Vector2d
r_hat (double xi, double eta, const Eigen::VectorXd &coord_vec, int poly_ord,
       std::string_view shape_function)
{
    Eigen::VectorXd unused;
    Eigen::MatrixXd shape_mat;
    std::tie (unused, shape_mat) = shape_N (shape_function, eta, poly_ord);
    return xi * shape_mat * coord_vec;
}

Eigen::Vector2d
r_c (double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
     const Eigen::Vector2d &centre, std::string_view shape_function_type)
{
    if (shape_function_type == "standard shape functions"
        || shape_function_type == "hierarchical shape functions")
        {
            return (r_hat_c (1, eta, poly_ord, coord_vec, centre) - centre);
        }
    else
        {
            std::cerr << "ERROR: No valid shape function type in FUNCTION "
                         "rc(eta,coord_vec,centre)"
                      << std::endl;
            return Eigen::Vector2d::Zero ();
        }
}

Eigen::Vector2d
r (double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
   const Eigen::Vector2d &centre, std::string_view shape_function_type)
{
    if (shape_function_type == "standard shape functions"
        || shape_function_type == "hierarchical shape functions")
        {
            return r_hat (1, eta, coord_vec, poly_ord);
        }
    else
        {
            std::cerr << "ERROR: No valid shape function type in FUNCTION "
                         "-r(eta,coord_vec)"
                      << std::endl;
            return Eigen::Vector2d::Zero ();
        }
}

Eigen::Matrix2d
j_mat (double eta, const Eigen::VectorXd &coord_vec, int poly_ord,
       std::string shape_function_type, const Eigen::Vector2d &centre)
{
    Eigen::Vector2d r_res
        = r (eta, poly_ord, coord_vec, shape_function_type, centre);
    Eigen::VectorXd shape_dN_res;
    Eigen::MatrixXd unused;
    std::tie (shape_dN_res, unused) = shape_dN (eta, poly_ord);
    Eigen::Vector2d dot_res = shape_dN_res[1] * coord_vec;
    return (Eigen::Matrix2d () << r_res[0], r_res[1], dot_res[0], dot_res[1])
        .finished ();
}
