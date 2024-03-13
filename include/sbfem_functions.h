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

/**
 * @brief Resize and pad a Eigen::VectorXd object with zeros.
 *
 * This function takes an Eigen::VectorXd object `vec` and an integer `n` as
 * input. It resizes the vector to a new size by adding `n` additional elements
 * at the end. The additional elements are set to zero.
 *
 * @param vec The original Eigen::VectorXd object.
 * @param n The number of elements to add at the end of the vector.
 * @return The resized Eigen::VectorXd object with zero-padded elements.
 */
Eigen::VectorXd resize_and_pad_with_zeros(const Eigen::VectorXd &vec, int n);

/**
 * @brief Calculates the position vector of a point on an element given its
 * local coordinates.
 *
 * This function takes the local coordinates (xi, eta) of a point on an element,
 * the order of the polynomial basis functions, the coordinates of the
 * Lagrangian points of the element, the type of shape function, and the centre
 * of the element as input. It calculates the position vector r_hat_c of the
 * point using Lagrange interpolation.
 *
 * @param xi              The xi coordinate of the point on the element.
 * @param eta             The eta coordinate of the point on the element.
 * @param poly_ord        The order of the polynomial basis functions.
 * @param coord_vec       The coordinates of the Lagrangian points of the
 * element.
 * @param shape_function  The type of shape function used for interpolation.
 * @param centre          The centre of the element.
 * @return                The position vector r_hat_c of the point.
 */
Eigen::Vector2d r_hat_c(
    double xi, double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
    ShapeFunctionType shape_function,
    const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

/**
 * @brief Calculates the r-hat vector for a given evaluation point on an
 * element.
 *
 * This function calculates the r-hat vector for a given evaluation point (xi,
 * eta) on an element. The r-hat vector is computed by multiplying the shape
 * functions (shape_N_res) by the coordinate vector (coord_vec_used), and then
 * scaling it by xi.
 *
 * @param xi           The xi coordinate of the evaluation point
 * @param eta          The eta coordinate of the evaluation point
 * @param poly_ord     The order of the polynomial basis functions used in the
 * shape function calculation
 * @param coord_vec    The coordinate vector
 * @param shape_function The type of shape function (STANDARD or HIERARCHICAL)
 * @return             The r-hat vector
 *
 * @see shape_N
 * @see resize_and_pad_with_zeros
 */
Eigen::Vector2d r_hat(double xi, double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function);

/**
 * @brief Calculates the position vector of a point on an element given its
 * shape function type and coordinates.
 *
 * This function calculates the position vector r_c of a point on an element
 * using Lagrange interpolation. The position vector is calculated based on the
 * provided shape function type, polynomial order, coordinates of Lagrangian
 * points, and the centre of the element.
 *
 * @param eta                  The local coordinate of the point on the element.
 * @param poly_ord             The order of the polynomial basis functions.
 * @param coord_vec            The coordinates of the Lagrangian points of the
 * element.
 * @param shape_function_type  The type of shape function used for
 * interpolation.
 * @param centre               The centre of the element.
 * @return                     The position vector r_c of the point.
 *
 * @throws                     std::cerr if the shape_function_type is not
 * valid.
 * @since                      1.0.0
 */
Eigen::Vector2d r_c(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                    ShapeFunctionType shape_function_type,
                    const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

/**
 * @brief Calculates the r vector for a given evaluation point on an element.
 *
 * This function calculates the r vector for a given evaluation point (eta) on
 * an element. The r vector is computed by calling the r_hat function with xi =
 * 1, which returns the r-hat vector. If the shape function type is STANDARD or
 * HIERARCHICAL, the r-hat vector is computed and returned. Otherwise, an error
 * message is printed and a zero vector is returned.
 *
 * @param eta             The eta coordinate of the evaluation point
 * @param poly_ord        The order of the polynomial basis functions used in
 * the shape function calculation
 * @param coord_vec       The coordinate vector
 * @param shape_function  The type of shape function (STANDARD or HIERARCHICAL)
 * @return                The r vector
 *
 * @see r_hat
 */
Eigen::Vector2d r(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                  ShapeFunctionType shape_function_type);

/**
 * @brief Calculates the j_hat matrix for a given xi, eta, polynomial order,
 * coordinates, shape function type, and centre.
 *
 * This function calculates the j_hat matrix based on the provided xi, eta,
 * polynomial order, coordinates of Lagrangian points, shape function type, and
 * centre of the element. The j_hat matrix is a 2x2 matrix that combines the
 * position vector and derivative values.
 *
 * @param xi                    The local coordinate of the point on the
 * element.
 * @param eta                   The local coordinate of the point on the
 * element.
 * @param poly_ord              The order of the polynomial basis functions.
 * @param coord_vec             The coordinates of the Lagrangian points of the
 * element.
 * @param shape_function_type   The type of shape function used for
 * interpolation.
 * @param centre                The centre of the element.
 *
 * @return                      The j_hat matrix.

 */
Eigen::Matrix2d j_hat_mat_c(
    double xi, double eta, int poly_ord, const Eigen::MatrixXd &coord_vec,
    ShapeFunctionType shape_function_type,
    const Eigen::VectorXd &centre = Eigen::VectorXd::Zero(2));

/**
 * @brief Calculates the J-hat matrix for a given evaluation point on an
 * element.
 *
 * This function calculates the J-hat matrix for a given evaluation point (xi,
 * eta) on an element. The J-hat matrix is computed based on the given xi and
 * eta values, polynomial order, coordinate vector, and shape function type. It
 * calls the r() function to calculate the r vector and the shape_dN() function
 * to calculate the shape function derivative matrix. If the shape function type
 * is HIERARCHICAL, the coordinate vector is resized and padded with zeros. The
 * J-hat matrix is created and returned as an Eigen::Matrix2d object.
 *
 * @param xi                   The xi coordinate of the evaluation point
 * @param eta                  The eta coordinate of the evaluation point
 * @param poly_ord             The order of the polynomial basis functions used
 * in the shape function calculation
 * @param coord_vec            The coordinate vector
 * @param shape_function_type  The type of shape function (STANDARD or
 * HIERARCHICAL)
 * @return                     The J-hat matrix
 *
 * @see r()
 * @see shape_dN()
 * @see resize_and_pad_with_zeros()
 */
Eigen::Matrix2d j_hat_mat(double xi, double eta, int poly_ord,
                          const Eigen::MatrixXd &coord_vec,
                          ShapeFunctionType shape_function_type);

/**
 * @brief Calculates the jacobian matrix for a given evaluation point on an
 * element.
 *
 * This function calculates the jacobian matrix for a given evaluation point
 * (eta) on an element. The function first computes the r vector using the r()
 * function. Then, it computes the shape function derivative matrix using the
 * shape_dN() function. Finally, depending on the shape function type, the
 * coordinate vector is resized and padded with zeros. The dot product of the
 * shape function derivative matrix and the coordinate vector is calculated and
 * stored in the dot_res vector. The jacobian matrix is then constructed using
 * the r_res vector (r() function result) and dot_res vector.
 *
 * @param eta               The eta coordinate of the evaluation point
 * @param poly_ord          The order of the polynomial basis functions used in
 * the shape function calculation
 * @param coord_vec         The coordinate vector
 * @param shape_function_type    The type of shape function (STANDARD or
 * HIERARCHICAL)
 * @return                  The jacobian matrix
 *
 * @see r, shape_dN, resize_and_pad_with_zeros
 */
Eigen::Matrix2d j_mat_c(
    double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
    ShapeFunctionType shape_function_type,
    const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

Eigen::Matrix2d j_mat(double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function_type);

/**
 * Calculates the determinant of the transformation matrix for a given value of
 * eta, polynomial order, coordinate vector, shape function type, and centre.
 *
 * @param eta The value of eta.
 * @param poly_ord The polynomial order.
 * @param coord_vec The coordinate vector.
 * @param shape_function_type The shape function type.
 * @param centre The centre point.
 * @return The determinant of the transformation matrix.
 */
double det_j(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
             ShapeFunctionType shape_function_type,
             const Eigen::Vector2d &centre = Eigen::Vector2d::Zero());

/**
 * @brief Calculates the derivative of the volumetric term (dV) yet to integrate
 * to xi and eta.
 *
 * This function calculates the derivative of the volumetric term (dV) yet to
 * integrate with respect to xi and eta, given the coordinates defined by xi and
 * eta, the polynomial order (poly_ord), the coordinate vector (coord_vec), and
 * the type of shape function (shape_function_type).
 *
 * @param xi The xi coordinate.
 * @param eta The eta coordinate.
 * @param poly_ord The polynomial order.
 * @param coord_vec The coordinate vector.
 * @param shape_function_type The type of shape function.
 * @return The derivative of the volumetric term (dV) yet to integrate with
 * respect to xi and eta.
 */
double dV_divided_by_dxi_deta(double xi, double eta, int poly_ord,
                              const Eigen::VectorXd &coord_vec,
                              ShapeFunctionType shape_function_type);


/**
 * @brief Calculates the transformed coordinates in xi space for a given eta
 * value.
 *
 * This function takes an eta value, polynomial order, coordinate vector, and
 * shape function type as input, and returns the transformed coordinates in xi
 * space as a column vector.
 *
 * @param eta The eta value for which to calculate the transformed coordinates.
 * @param poly_ord The polynomial order.
 * @param coord_vec The coordinate vector.
 * @param shape_function_type The shape function type.
 * @return The transformed coordinates in xi space.
 */
Eigen::Vector2d g_xi(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                     ShapeFunctionType shape_function_type);

/**
 * @brief Calculates the unit normal vector in xi space for a given eta value.
 *
 * This function takes an eta value, polynomial order, coordinate vector, and
 * shape function type as input, and returns the unit normal vector in xi space
 * as a 2-dimensional vector.
 *
 * @param eta The eta value for which to calculate the unit normal vector.
 * @param poly_ord The polynomial order.
 * @param coord_vec The coordinate vector.
 * @param shape_function_type The shape function type.
 * @return The unit normal vector in xi space.
 */
Eigen::Vector2d n_xi(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                     ShapeFunctionType shape_function_type);


/**
 * @brief Calculate the g_eta function.
 *
 * This function calculates the g_eta function, which is used in certain
 * mathematical calculations. It takes the input parameters eta, poly_ord,
 * coord_vec, and shape_function_type, and returns a 2D Eigen::Vector2d
 * containing the result. The function utilizes the r function to calculate an
 * intermediate result, which is then used to calculate the final result.
 *
 * @param eta The value of eta parameter.
 * @param poly_ord The value of poly_ord parameter.
 * @param coord_vec The input vector of coordinates.
 * @param shape_function_type The type of shape function to be used.
 * @return Eigen::Vector2d The resulting 2D vector.
 *
 * @see ShapeFunctionType
 * @see r
 */
Eigen::Vector2d g_eta(double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function_type);


/**
 * @brief Calculate the n_eta function.
 *
 * This function calculates the n_eta function, which is used in certain
 * mathematical calculations. It takes the input parameters eta, poly_ord,
 * coord_vec, and shape_function_type, and returns a 2D Eigen::Vector2d
 * containing the result. The function utilizes the g_eta function to calculate
 * an intermediate result, which is then normalized to obtain the final result.
 *
 * @param eta The value of the eta parameter.
 * @param poly_ord The value of the poly_ord parameter.
 * @param coord_vec The input vector of coordinates.
 * @param shape_function_type The type of shape function to be used.
 * @return Eigen::Vector2d The resulting 2D vector.
 *
 * @see Eigen::Vector2d
 * @see ShapeFunctionType
 * @see g_eta
 */
Eigen::Vector2d n_eta(double eta, int poly_ord,
                      const Eigen::VectorXd &coord_vec,
                      ShapeFunctionType shape_function_type);



Eigen::MatrixXd b1(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre);


Eigen::MatrixXd b2(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre);

Eigen::MatrixXd B1(double eta, int poly_ord, const Eigen::VectorXd &coord_vec,
                   ShapeFunctionType shape_function_type,
                   const Eigen::Vector2d &centre);

#endif // SBFEM_SBFEM_FUNCTIONS_H
