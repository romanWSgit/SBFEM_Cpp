//
// Created by Roman Wallner- Silberhuber on 06.03.24.
//

#ifndef FORTRAN_INTEROPERABILITY_H
#define FORTRAN_INTEROPERABILITY_H

#include <memory>
#include <vector>

// Declare Fortran functions for use in C++
extern "C" {
    void n_vec_f_c(double eta, int poly_ord, bool deriv, double* n_vec);
}


// Declaration of utility function (if used in multiple places or for clarity)
std::vector<double> callNVecF(double eta, int poly_ord, bool deriv);

#endif // FORTRAN_INTEROPERABILITY_H
