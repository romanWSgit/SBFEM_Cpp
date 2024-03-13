//
// Created by Roman Wallner- Silberhuber on 06.03.24.
//
#include "FortranInteroperability.h"
#include <iostream>



// Implement the utility function
// std::vector<double> callNVecF(double eta, int polyOrd, bool deriv)
//{
//    int n;
//    auto rawPtr = static_cast<double *>(n_vec_f_c(eta, polyOrd, deriv, n));
//    std::unique_ptr<double[], FortranArrayDeleter> smartPtr(
//        rawPtr, FortranArrayDeleter(n));
//    return std::vector<double>(
//        smartPtr.get(), smartPtr.get() + n); // Convert to std::vector<double>
//}

std::vector<double> callNVecF(double eta, int poly_ord, bool deriv)
{
    std::vector<double> n_vec(poly_ord +
                              1); // Allocate vector with size poly_ord + 1
    n_vec_f_c(eta, poly_ord, deriv,
              n_vec.data()); // Pass the underlying array to Fortran

    return n_vec; // n_vec is now filled by the Fortran subroutine
}
