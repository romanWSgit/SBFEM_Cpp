//
// Created by Roman Wallner- Silberhuber on 04.06.23.
//

#ifndef SBFEM_SBFEM_PLOT_H
#define SBFEM_SBFEM_PLOT_H

#include "SuperElementJson.h"
#include "gnuplot-iostream.h"
#include <Eigen/Dense>
#include <matplot/matplot.h>

void plotSbfemSuperelement(const SuperElementJson & sE, bool gnuplot=true);

void visualizeDoubleMatrix (const Eigen::MatrixXd &matrix);

void visualizeOnlyDiagonalOfDoubleMatrix (const Eigen::MatrixXd &matrix);

#endif //SBFEM_SBFEM_PLOT_H
