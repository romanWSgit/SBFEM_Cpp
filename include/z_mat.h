//
// Created by Roman Wallner- Silberhuber on 26.01.24.
//

#ifndef SBFEM_Z_MAT_H
#define SBFEM_Z_MAT_H

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>



/**
 * @class NormE0
 * @brief Class representing a norm E0.
 *
 * The NormE0 class represents a norm E0, which consists of an index and a
 * factor matrix.
 */
struct NormE0Data
{
    Eigen::Vector2i idx;
    Eigen::MatrixXd fctr;
};


/**
 * @brief Computes the norm E0 of a matrix.
 *
 * This function computes the norm E0 of a given matrix. The norm E0 is
 * calculated for each node in the matrix and returns a vector of NormE0Data
 * objects, where each object represents the index and factor matrix of a
 * specific node.
 *
 * @param E0 The input matrix.
 * @param nnodes The number of nodes in the matrix.
 * @param ndn The number of degrees of freedom per node.
 * @return A vector of NormE0Data objects representing the norm E0 of each
 * node.
 */
std::vector<NormE0Data> computeNormE0 (const Eigen::MatrixXd &E0, int nnodes,
                                       int ndn);

/**
 * @brief Normalizes the data in the NormE0Array using the values in E0, E1,
 * and E2.
 *
 * The function adjusts the values in the NormE0Array based on the provided E0,
 * E1, and E2 vectors. Each element in the NormE0Array is divided by the
 * corresponding elements in the E0, E1, and E2 vectors, resulting in
 * normalized values.
 *
 * @param NormE0Array The array of NormE0Data elements to be normalized.
 * @param E0 The vector of values for E0.
 * @param E1 The vector of values for E1.
 * @param E2 The vector of values for E2.
 */
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
normalizeE0E1E2 (std::vector<NormE0Data> &NormE0Array, Eigen::MatrixXd &E0,
                 Eigen::MatrixXd &E1, Eigen::MatrixXd &E2);

/**
 * @brief Calculate the inverse norm of E0M0.
 *
 * This function calculates the inverse norm of E0M0 given an array of
 * NormE0Data and a vector v. The calculation is based on the following
 * formula:
 *
 *                  nd
 *          --------------------
 *         \| (v[i] - NormE0Array[i].E0) ^ 2
 * i = 0
 *
 * @param NormE0Array The array of NormE0Data.
 * @param v The vector v.
 * @param nd The size of the NormE0Array and v.
 * @return The inverse norm of E0M0.
 */
Eigen::MatrixXd computeInverseNormE0M0 (std::vector<NormE0Data> &NormE0Array,
                                  Eigen::MatrixXd &v, int nd);

/**
 * @brief Calculates the condition number of a given matrix.
 * 
 * The condition number is a measure of how sensitive the output values are
 * to changes in the input values. It is defined as the ratio of the largest
 * singular value to the smallest singular value of the matrix. A higher
 * condition number indicates a higher sensitivity to input changes.
 * 
 * @param matrix A const reference to the matrix for which to calculate the condition number.
 * 
 * @return The condition number of the matrix.
 */
double calculateConditionNumber(const Eigen::MatrixXd &matrix);


/**
 * \fn Zmatrix(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2,
 * int ndn) \brief Constructs a Z-matrix using the given input matrices and
 * dimensions.
 *
 * This function constructs a Z-matrix using the provided input matrices and
 * dimensions.
 *
 * \param E0 The first input matrix.
 * \param E1 The second input matrix.
 * \param E2 The third input matrix.
 * \param ndn The number of degrees of freedom.
 */
Eigen::MatrixXd matrixZMethod (Eigen::MatrixXd &E0, Eigen::MatrixXd &E1,
                         Eigen::MatrixXd &E2, int ndn);

#endif // SBFEM_Z_MAT_H
