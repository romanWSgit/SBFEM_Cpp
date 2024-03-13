/**
 * \file   sbfem_functions.cpp
 * \author Roman Wallner- Silberhuber
 * \date   04.06.2023
 * \brief  This file contains the SBFEM - transformation of geometry functions
 */

#include "sbfem_functions.h"
#include "sbfem_math.h"



Eigen::VectorXd resize_and_pad_with_zeros(const Eigen::VectorXd &vec, int n)
{
    Eigen::VectorXd resized = vec;
    resized.conservativeResize(vec.size() + n);
    for (int i = 1; i <= n; ++i)
    {
        resized[resized.size() - i] = 0;
    }
    return resized;
}

Eigen::Vector2d r_hat_c(double xi, double eta, int poly_ord,
                        const Eigen::VectorXd &coord_vec,
                        ShapeFunctionType shape_function,
                        const Eigen::Vector2d &centre)
{
    Eigen::MatrixXd shape_N_res;
    std::tie(std::ignore, shape_N_res) = shape_N(eta, poly_ord, shape_function);

    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }

    return xi * shape_N_res * coord_vec_used + centre;
}

Eigen::Vector2d r_hat(double xi, double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function)
{

    Eigen::MatrixXd shape_N_res;
    std::tie(std::ignore, shape_N_res) = shape_N(eta, poly_ord, shape_function);

    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }

    return xi * shape_N_res * coord_vec_used;
}

Eigen::Vector2d r_c(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                    ShapeFunctionType shape_function_type,
                    const Eigen::Vector2d &centre)
{
    if (shape_function_type == ShapeFunctionType::STANDARD ||
        shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        return (
            r_hat_c(1, eta, poly_ord, coord_vec, shape_function_type, centre) -
            centre);
    }
    else
    {
        std::cerr << "ERROR: No valid shape function type in FUNCTION "
                     "rc(eta,coord_vec,centre)"
                  << std::endl;
        return Eigen::Vector2d::Zero();
    }
}

Eigen::Vector2d r(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                  ShapeFunctionType shape_function_type)
{
    if (shape_function_type == ShapeFunctionType::STANDARD ||
        shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        return r_hat(1, eta, poly_ord, coord_vec, shape_function_type);
    }
    else
    {
        std::cerr << "ERROR: No valid shape function type in FUNCTION "
                     "-r(eta,coord_vec)"
                  << std::endl;
        return Eigen::Vector2d::Zero();
    }
}

Eigen::Matrix2d j_hat_mat_c(double xi, double eta, int poly_ord,
                            const Eigen::MatrixXd &coord_vec,
                            ShapeFunctionType shape_function_type,
                            const Eigen::VectorXd &centre)
{
    Eigen::Vector2d r_res =
        r_c(eta, poly_ord, coord_vec, shape_function_type, centre);

    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }

    Eigen::MatrixXd result(2, 2);

    Eigen::Vector2d dot_res = xi * shape_dN_res * coord_vec_used;
    result(0, 0) = r_res(0);
    result(0, 1) = r_res(1);
    result(1, 0) = dot_res(0);
    result(1, 1) = dot_res(1);

    return result;
}

Eigen::Matrix2d j_hat_mat(double xi, double eta, int poly_ord,
                          const Eigen::MatrixXd &coord_vec,
                          ShapeFunctionType shape_function_type)
{
    Eigen::Vector2d r_res = r(eta, poly_ord, coord_vec, shape_function_type);

    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }

    Eigen::MatrixXd result(2, 2);

    Eigen::Vector2d dot_res = xi * shape_dN_res * coord_vec_used;
    result(0, 0) = r_res(0);
    result(0, 1) = r_res(1);
    result(1, 0) = dot_res(0);
    result(1, 1) = dot_res(1);

    return result;
}

Eigen::Matrix2d j_mat_c(double eta, int poly_ord,
                        const Eigen::VectorXd &coord_vec,
                        ShapeFunctionType shape_function_type,
                        const Eigen::Vector2d &centre)
{
    Eigen::Vector2d r_res =
        r_c(eta, poly_ord, coord_vec, shape_function_type, centre);

    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }
    Eigen::Vector2d dot_res = shape_dN_res * coord_vec_used;
    return (Eigen::Matrix2d() << r_res[0], r_res[1], dot_res[0], dot_res[1])
        .finished();
}

Eigen::Matrix2d j_mat(double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function_type)
{
    Eigen::Vector2d r_res = r(eta, poly_ord, coord_vec, shape_function_type);

    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }
    Eigen::Vector2d dot_res = shape_dN_res * coord_vec_used;
    return (Eigen::Matrix2d() << r_res[0], r_res[1], dot_res[0], dot_res[1])
        .finished();
}

double det_j(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
             ShapeFunctionType shape_function_type,
             const Eigen::Vector2d &centre)
{

    double det =
        j_mat(eta, poly_ord, coord_vec, shape_function_type).determinant();
    return det;
}

double dV_divided_by_dxi_deta(double xi, double eta, int poly_ord,
                              const Eigen::VectorXd &coord_vec,
                              ShapeFunctionType shape_function_type)
{
    return xi * det_j(eta, poly_ord, coord_vec, shape_function_type);
}

// normal Vectors:

// g\[Xi][\[Eta]_] = -{-D[r[\[Eta]][[2]], \[Eta]],
//                      D[r[\[Eta]][[1]], \[Eta]]} // FullSimplify;
//                    g\[Eta][\[Eta]_] = {-r[\[Eta]][[2]], r[\[Eta]][[1]]} //
//                    Simplify;
// n\[Xi][\[Eta]_] =
//     g\[Xi][\[Eta]]/Sqrt[g\[Xi][\[Eta]] . g\[Xi][\[Eta]]] // Simplify;
//     n\[Eta][\[Eta]_] =
//         g\[Eta][\[Eta]]/Sqrt[g\[Eta][\[Eta]] . g\[Eta][\[Eta]]] // Simplify;

// Example for g_xi, similar adaptations needed for g_eta, n_xi, n_eta, b1, b2

Eigen::Vector2d g_xi(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                     ShapeFunctionType shape_function_type)
{

    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    Eigen::VectorXd int_result;
    Eigen::Vector2d result;

    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }
    int_result = shape_dN_res * coord_vec_used;
    result(0) = -int_result(1);
    result(1) = int_result(0);
    return result;
}

Eigen::Vector2d n_xi(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                     ShapeFunctionType shape_function_type)
{
    Eigen::Vector2d g = g_xi(eta, poly_ord, coord_vec, shape_function_type);
    Eigen::Vector2d n = g.normalized();
    return n;
}

Eigen::Vector2d g_eta(double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function_type)
{
    Eigen::Vector2d result;
    Eigen::Vector2d int_result =
        r(eta, poly_ord, coord_vec, shape_function_type);
    result(0) = -int_result(1);
    result(1) = int_result(0);
    return result;
}

Eigen::Vector2d n_eta(double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function_type)
{
    Eigen::Vector2d g = g_eta(eta, poly_ord, coord_vec, shape_function_type);
    Eigen::Vector2d n = g.normalized();
    return n;
}

// Todo: Implement b1, b2 and Dr
Eigen::MatrixXd b1(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre)
{
    double determinant =
        det_j(eta, poly_ord, coord_vec, shape_function_type, centre);
    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    Eigen::VectorXd int_result;

    Eigen::VectorXd coord_vec_used = coord_vec;
    if (shape_function_type == ShapeFunctionType::HIERARCHICAL)
    {
        coord_vec_used =
            resize_and_pad_with_zeros(coord_vec, (poly_ord - 1) * 2);
    }
    int_result = shape_dN_res * coord_vec_used;
    Eigen::MatrixXd result(3, 2);

    result << int_result(1) / determinant, 0, 0, -int_result(0) / determinant,
        -int_result(0) / determinant, int_result(1) / determinant;

    return result;
}

Eigen::MatrixXd b2(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre)
{
    double determinant =
        det_j(eta, poly_ord, coord_vec, shape_function_type, centre);
    Eigen::Vector2d int_result =
        r(eta, poly_ord, coord_vec, shape_function_type);
    Eigen::MatrixXd result(3, 2);

    result << -int_result(1) / determinant, 0, 0, int_result(0) / determinant,
        int_result(0) / determinant, -int_result(1) / determinant;

    return result;
}

Eigen::MatrixXd B1(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre)
{

    Eigen::MatrixXd b1mat =
        b1(eta, poly_ord, coord_vec, shape_function_type, centre);
    Eigen::MatrixXd shape_N_res;
    std::tie(std::ignore, shape_N_res) =
        shape_N(eta, poly_ord, shape_function_type);
    // bi has the shape 3x2 and N the shape
    return b1mat * shape_N_res;
}


Eigen::MatrixXd B2(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre)
{

    Eigen::MatrixXd b2mat =
        b2(eta, poly_ord, coord_vec, shape_function_type, centre);
    Eigen::MatrixXd shape_dN_res;
    std::tie(std::ignore, shape_dN_res) =
        shape_dN(eta, poly_ord, shape_function_type);
    // bi has the shape 3x2 and N the shape
    return b2mat * shape_dN_res;
}

