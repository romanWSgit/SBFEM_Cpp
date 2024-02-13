/**
 * \file   helper_functions.h
 * \author Roman Wallner- Silberhuber
 * \date   25.05.2023
 * \brief  This head contains some usefull helper functions
 */

#ifndef SBFEM_HELPER_FUNCTIONS_H
#define SBFEM_HELPER_FUNCTIONS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <nlohmann/json.hpp>
#include <tuple>
#include <algorithm>
#include <optional>
#include <vector>
#include <fstream>
#include "MaterialData.h"
#include <iostream>
#include <type_traits>
#include "fast_matrix_market/fast_matrix_market.hpp"
#include "fast_matrix_market/app/Eigen.hpp"

/**

   * @brief Writes an Eigen matrix to a Matrix Market (.mtx) file.

   * The Matrix Market format is widely used for representing matrices in numerical computations.
     !!! NOTE Template needs to be in the header
   * @param matrix The Eigen matrix to be written.

   * @param filename The name of the .mtx file to be generated.

   * @details This function takes an Eigen matrix and writes it to a Matrix Market (.mtx) file.

   * The Matrix Market format is a widely-used file format for representing matrices in numerical computations.

   * The generated file will have a header specifying the matrix details, followed by the matrix data.

   * @throws std::ofstream::failure If the file cannot be opened.

   * @return void

   */
template<class Matrix>
void writeEigenMatrixToMtx (const Matrix &matrix, const std::string &filename)
{
    std::ofstream file (filename);
    if (file.is_open ())
        {
            file << "%%MatrixMarket matrix array real general\n";
            file << matrix.rows () << " " << matrix.cols () << "\n";
            for (int i = 0; i < matrix.rows (); ++i)
                {
                    for (int j = 0; j < matrix.cols (); ++j)
                        {
                            file << matrix (i, j) << "\n";
                        }
                }
            file.close ();
        }
    else
        {
            std::cout << "Unable to open file";
        }
}

/**
 * @brief Reads a Matrix Market file into an Eigen Dense matrix.
 *
 * This function takes a file path as input and reads a Matrix Market file into an Eigen Dense matrix.
 * The function uses the fast_matrix_market library to handle the Matrix Market format.
 *
 * @param filePath The file path of the Matrix Market file to be read.
 * @return The Eigen Dense matrix containing the data from the Matrix Market file.
 * @throw invalid_mm If the Matrix Market file is invalid or cannot be read.
 */
Eigen::MatrixXd readMatrixFromFile (const std::string &filePath);

/**
 * @brief Function to convert a JSON array to an Eigen matrix.
 *
 * This function takes a JSON array as input and converts it into an Eigen matrix. The JSON array is assumed to represent a 2D matrix, where each element in the array is a row of the
* matrix.
 *
 * @param jsonArray The JSON array to be converted.
 * @return The Eigen matrix created from the JSON array.
 *
 * @note The JSON array should be a valid 2D matrix representation.
 * @warning The function assumes that all rows in the JSON array have the same number of elements.
 *
 * @see fromJsonArray()
 */
template<typename MatrixType>
MatrixType jsonToEigenMatrix (const nlohmann::json &jsonArray)
{
    // Check if the input is a non-empty array of arrays
    if (jsonArray.empty () || !jsonArray.is_array ()
        || !jsonArray[0].is_array ())
        {
            throw std::invalid_argument ("Invalid JSON input: Expected a non-empty array of arrays.");
        }

    size_t rows = jsonArray.size ();
    size_t cols = jsonArray[0].size ();

    // Check for consistency in each row
    for (const auto &row: jsonArray)
        {
            if (!row.is_array () || row.size () != cols)
                {
                    throw std::invalid_argument ("Invalid JSON input: Inconsistent array sizes in matrix.");
                }
        }

    // Instantiate an Eigen matrix of the appropriate type
    MatrixType mat (rows, cols);
    for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < cols; ++j)
                {
                    // Use std::is_same to check the matrix type and cast accordingly
                    if constexpr (std::is_same<MatrixType, Eigen::MatrixXd>::value)
                        {
                            mat (i, j) = jsonArray[i][j].get<double> ();
                        }
                    else if constexpr (std::is_same<MatrixType, Eigen::MatrixXi>::value)
                        {
                            mat (i, j) = jsonArray[i][j].get<int> ();
                        }
                }
        }
    return mat;
}


/**
 * @brief Convert an Eigen vector to a std::vector
 *
 * This function takes an Eigen vector object and converts it to a std::vector<double>.
 * The resulting std::vector will have the same size as the input Eigen vector.
 *
 * @tparam Derived The derived Eigen vector type
 * @param vector The input Eigen vector to be converted
 * @return The resulting std::vector<double> object
 *
 * @note The input Eigen vector type must be an Eigen::MatrixBase derived class.
 */
template<typename Derived>
std::vector<double> eigenToStdVector (const Eigen::MatrixBase<Derived> &vector)
{
    std::vector<double> vec (vector.size ());
    for (int i = 0; i < vector.size (); ++i)
        {
            vec[i] = static_cast<double>(vector (i));
        }
    return vec;
}


/**
* \brief Convert an Eigen matrix to a 2D std::vector
*
* This function converts an Eigen matrix to a 2D std::vector with double precision.
* The resulting 2D std::vector has the same dimensions as the input matrix.
*
* \tparam Derived The template parameter type for the Eigen matrix.
* \param matrix The Eigen matrix to be converted.
* \return A 2D std::vector with the same dimensions as the input matrix.
*/
template<typename Derived>
std::vector<std::vector<double>>
eigenToStdVector2D (const Eigen::MatrixBase<Derived> &matrix)
{
    std::vector<std::vector<double>> vec (matrix.rows (), std::vector<double> (matrix.cols ()));
    for (int i = 0; i < matrix.rows (); ++i)
        {
            for (int j = 0; j < matrix.cols (); ++j)
                {
                    vec[i][j] = static_cast<double>(matrix (i, j));
                }
        }
    return vec;
}

/**
 * \fn     StructData get_material_data_for_type(std::string& material_name, \
                                      std::vector<StructData>& material_data_list);
 * \brief  Short description of the function
 * \param  nlohmann::json& j Description of the first parameter
 * \param  std::string field Description of the second parameter
 * \return returns an Eigen::MatrixXd matrix iff the field name is found in the json handle
 *
 * More detailed description of the function
*/
StructData getMaterialDataForType (std::string &material_name, \
                                      std::vector<StructData> &material_data_list);

#endif //SBFEM_HELPER_FUNCTIONS_H
