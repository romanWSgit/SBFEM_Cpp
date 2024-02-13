//
// Created by Roman Wallner- Silberhuber on 13.07.23.
//

#ifndef SBFEM_REORDER_SCHUR_H
#define SBFEM_REORDER_SCHUR_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/LU>
#include <vector>
#include <utility>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>

struct SchurData
{
    Eigen::MatrixXd U;
    Eigen::MatrixXd S;
    Eigen::VectorXd diagonalElements;
    double maxDiagonalElement;
    double minDiagonalElement;
    Eigen::VectorXcd eigenvaluesOfS;
    std::complex<double> maxEV;
    std::complex<double> minEV;

    explicit SchurData(const Eigen::RealSchur<Eigen::MatrixXd>& A);

    SchurData(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& S);

    ~SchurData();

    void getMinMaxEv();


    SchurData normalizeInPlace(SchurData& schurUS, const Eigen::VectorXi& v);

};

std::vector<int> swaplist(Eigen::VectorXcd p, std::vector<int> s, std::complex<double> z, double b);

std::pair<double, Eigen::Index> select(const Eigen::VectorXcd& p, std::complex<double> z);

SchurData normalize(SchurData schurUS, Eigen::VectorXi v);

Eigen::Matrix2d rot(const Eigen::Matrix2d& X);

void swap(SchurData schurData, Eigen::VectorXi v, Eigen::VectorXi w);

Eigen::VectorXd SRSchur(SchurData& schurData, double z, int b);

#endif //SBFEM_REORDER_SCHUR_H
