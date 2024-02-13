/**
 * \file    main.cpp
 * \author  Roman Wallner- Silberhuber
 * \date    21.05.23
 * \brief   The main method moderates the Sbfem Computation
 *
 * \details One thing to note here...
 * \todo    add details to the File discripton
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/version.hpp>
#include <matplot/matplot.h> // For Matplot++

#include "MaterialData.h"
#include "SuperElement.h"
#include "sbfem_math.h"

#include "helper_functions.h"
#include "plot.h"
#include "sbfem_functions.h"

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <nlohmann/json.hpp>
#include <utility>
#include <vector>

#include "fast_matrix_market/app/Eigen.hpp"
#include "fast_matrix_market/fast_matrix_market.hpp"

// Own
#include "reorder_schur.h"
#include "z_mat.h"

extern "C"
{
    void
    mb03qd_ (char *DICO, char *STDOM, char *JOBU, int *N, int *NLOW, int *NSUP,
             double *ALPHA, double *A, int *LDA, double *U, int *LDU, int *NDIM,
             double *DWORK, int *INFO);
}

constexpr bool IMPORT_SCHUR_FLAG = true;

int
main ()
{
    Eigen::MatrixXd matrixZ;
    Eigen::MatrixXd matrixQPython;
    Eigen::MatrixXd matrixRPython;

    Eigen::MatrixXd matrixZImported = readMatrixFromFile (
        "/Users/roman_w_s/Developer/PYTHON/SBFEM/schur/data/Zmat.mtx");
    Eigen::MatrixXd matrixE0Imported = readMatrixFromFile (
        "/Users/roman_w_s/Developer/PYTHON/SBFEM/schur/data/E0.mtx");
    Eigen::MatrixXd matrixE1Imported = readMatrixFromFile (
        "/Users/roman_w_s/Developer/PYTHON/SBFEM/schur/data/E1.mtx");
    Eigen::MatrixXd matrixE2Imported = readMatrixFromFile (
        "/Users/roman_w_s/Developer/PYTHON/SBFEM/schur/data/E2.mtx");

    int numberOfRows = matrixE0Imported.rows ();

    std::vector normEO = computeNormE0 (matrixE0Imported, numberOfRows / 2, 2);

    auto [normalizedE0, normalizedE1, normalizedE2] = normalizeE0E1E2 (
        normEO, matrixE0Imported, matrixE1Imported, matrixE2Imported);

    Eigen::MatrixXd matrixZNormed
        = matrixZMethod (normalizedE0, normalizedE1, normalizedE2, 2);
    Eigen::RealSchur<Eigen::MatrixXd> schurZ (matrixZImported);
    Eigen::RealSchur<Eigen::MatrixXd> schurZNormed (matrixZNormed);

    Eigen::MatrixXd matrixUSchurNormed = schurZNormed.matrixU ();
    Eigen::MatrixXd matrixTSchurNormed = schurZNormed.matrixT ();
    Eigen::MatrixXd matrixUSchur = schurZ.matrixU ();
    Eigen::MatrixXd matrixTSchur = schurZ.matrixT ();
    std::cout << "U * T * U^T = " << std::endl
              << matrixUSchurNormed * matrixTSchurNormed
                         * matrixUSchurNormed.transpose ()
                     - matrixZNormed
              << std::endl;

    double conditionNumberNormed
        = calculateConditionNumber (matrixUSchurNormed);
    std::cout << "Condition Number of normalized U-matrix: "
              << conditionNumberNormed << std::endl;

    double conditionNumber = calculateConditionNumber (matrixUSchur);
    std::cout << "Condition Number of U-matrix: " << conditionNumber
              << std::endl;
    // visualizeDoubleMatrix (matrixTSchurNormed);
    // visualizeOnlyDiagonalOfDoubleMatrix (matrixTSchurNormed);

    std::unique_ptr<SchurData> schurDataPointer;
    if (IMPORT_SCHUR_FLAG)
        {
            std::ifstream file1 ("/Users/roman_w_s/Desktop/Q_schur.mtx");
            fast_matrix_market::read_matrix_market_eigen_dense (file1,
                                                                matrixQPython);
            file1.close ();
            std::ifstream file2 ("/Users/roman_w_s/Desktop/R_schur.mtx");
            fast_matrix_market::read_matrix_market_eigen_dense (file2,
                                                                matrixRPython);
            file2.close ();
            // Use std::make_unique to create a new SchurData object
            schurDataPointer
                = std::make_unique<SchurData> (matrixQPython, matrixRPython);
            schurDataPointer->getMinMaxEv ();
        }
    else
        {
            // Compute the Schur decomposition
            Eigen::RealSchur<Eigen::MatrixXd> A (matrixZImported);
            schurDataPointer = std::make_unique<SchurData> (A);
            schurDataPointer->getMinMaxEv ();
        }

    //    Eigen::VectorXd ap
    //        = SRSchur (*sE_Schur, sE_Schur->minDiagonalElement - 1.0, 0);

    //    // The matrix U in the decomposition is orthogonal
    //    std::cout << "The orthogonal matrix U in the decomposition of m
    //    is:\n" << U << std::endl; Eigen::Matrix<double, -1, 1> ap;
    //    std::vector<int> swapEXP;

    // Read JSON files
    std::ifstream domain_file (
        "/Users/roman_w_s/Desktop/SBFEM_DATA/rect_1.json");
    std::ifstream material_file ("/Users/roman_w_s/Developer/WL/Applications/"
                                 "SbfemPkg/Resources/Database/ds.json");
    nlohmann::json j_geom;
    nlohmann::json j_material;
    domain_file >> j_geom;
    material_file >> j_material;

    auto sE = SuperElement (j_geom);

    // Print Eigen array
    std::cout << sE.getMSEElementsM () << std::endl;
    std::cout << sE.getMSEShapeFct () << std::endl;

    // todo: write global try catch block
    try
        {
            auto [shapeVec, shapeMat]
                = shape_N (sE.getMSEShapeFct (), 0.1, sE.getMSEPolyOrd ());
            std::cout << shapeVec << '\n';
            std::cout << shapeMat << '\n';
        }
    catch (const std::runtime_error &e)
        {
            std::cout << e.what () << '\n';
        }

    std::vector<StructData> material_data_list;

    for (const auto &struct_json : j_material)
        {
            material_data_list.push_back (to_struct_data (struct_json));
        }

    // Now you can find a struct by its type:
    std::string search_materialType = "Flat Disc";

    StructData material
        = getMaterialDataForType (search_materialType, material_data_list);

    auto material_data
        = findStructDataAndMaterial (material_data_list, "Flat Disc", "steel");

    std::cout << material_data->params.Em << std::endl;

    std::cout << "ltg: \n" << sE.getMSELtg () << std::endl;

    plotSbfemSuperelement (sE, false);
    plotSbfemSuperelement (sE);

    std::cout << "Hello SBFEM!" << std::endl;

    std::cout << "Using Boost " << BOOST_VERSION / 100000
              << "."                               // major version
              << BOOST_VERSION / 100 % 1000 << "." // minor version
              << BOOST_VERSION % 100               // patch level
              << std::endl;

    return 0;
}
