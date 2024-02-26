//
// Created by Roman Wallner- Silberhuber on 26.02.24.
//
#include "sbfem_math.h"
#include "sbfem_driver.h"
#include "sbfem_functions.h"
#include <Eigen/Dense>
#include <omp.h> // opem mp
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <nlohmann/json.hpp>
#include "SuperElementJson.h"
#include <fstream>
#include <iostream>



BOOST_AUTO_TEST_SUITE(SBFEM_MP_AREA_TEST)

BOOST_AUTO_TEST_CASE(test_mp_area)
{

    // Read JSON files
    std::ifstream domain_file(
        "/Users/roman_w_s/Desktop/SBFEM_DATA/rect_1.json");

    nlohmann::json j_geom;

    domain_file >> j_geom;


    auto sE = SuperElementJson(j_geom);

    auto start = std::chrono::high_resolution_clock::now();

    omp_set_num_threads(10);


    double area{0};
#pragma omp parallel for reduction(                                            \
        + : area) // This line enables parallel for loop
    for (int i = 0; i < sE.getMSEIelem(); ++i)
    {


        Eigen::Matrix<double, 4, 1> elem = sE.getMSEElementsM().row(i);
        auto f_bind = std::bind(dV_divided_by_dxi_deta, std::placeholders::_1,
                                std::placeholders::_2, 1, elem,
                                ShapeFunctionType::STANDARD);

        area += integrateScalarKronrod2D(f_bind, 0.0, 1.0, -1.0, 1.0);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "area: " << area << "\n";
    std::cout << "Time taken by area calculation: " << duration.count()
              << " microseconds" << std::endl;
    std::cout << "Time taken by area calculation : "
              << ((double)duration.count()) / 1000000 << " seconds"
              << std::endl;
    double expected_output = 144;
    BOOST_CHECK_CLOSE(area, expected_output, 0.001);
}

BOOST_AUTO_TEST_SUITE_END()
