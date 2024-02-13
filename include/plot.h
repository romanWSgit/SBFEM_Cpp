//
// Created by Roman Wallner- Silberhuber on 04.06.23.
//

#ifndef SBFEM_SBFEM_PLOT_H
#define SBFEM_SBFEM_PLOT_H

#include <Eigen/Dense>
#include "gnuplot-iostream.h"
#include "SuperElement.h"
#include <matplot/matplot.h>


void plotSbfemSuperelement(const SuperElement& sE, bool gnuplot=true);

void visualizeDoubleMatrix (const Eigen::MatrixXd &matrix);

void visualizeOnlyDiagonalOfDoubleMatrix (const Eigen::MatrixXd &matrix);

#endif //SBFEM_SBFEM_PLOT_H
