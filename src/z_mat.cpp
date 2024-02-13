//
// Created by Roman Wallner- Silberhuber on 26.01.24.
//
#include "z_mat.h"

std::vector<NormE0Data>
computeNormE0 (const Eigen::MatrixXd &E0, int nnodes, int ndn)
{
    std::vector<NormE0Data> NormE0Array;

    for (int ii = 1; ii <= nnodes; ++ii)
        {
            int i1 = (ii - 1) * ndn;
            int i2 = ii * ndn - 1;

            Eigen::MatrixXd matrixBlock
                = E0.block (i1, i1, i2 - i1 + 1, i2 - i1 + 1);
            Eigen::EigenSolver<Eigen::MatrixXd> solver (matrixBlock);
            Eigen::MatrixXd vct = solver.eigenvectors ().real ();
            Eigen::MatrixXd d = vct.transpose () * matrixBlock * vct;

            Eigen::MatrixXd diag = d.diagonal ()
                                       .array ()
                                       .inverse ()
                                       .sqrt ()
                                       .matrix ()
                                       .asDiagonal ();
            Eigen::MatrixXd fctr = vct * diag;

            NormE0Array.push_back (
                NormE0Data{ Eigen::Vector2i (i1, i2), fctr });
        }

    return NormE0Array;
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
normalizeE0E1E2 (std::vector<NormE0Data> &NormE0Array, Eigen::MatrixXd &E0,
                 Eigen::MatrixXd &E1, Eigen::MatrixXd &E2)
{
    for (auto &ii1 : NormE0Array)
        {
            int i1 = ii1.idx (0);
            int i2 = ii1.idx (1);
            Eigen::MatrixXd fctr = ii1.fctr;

            E0.block (0, i1, E0.rows (),
                      i2 - i1 + 1) *= fctr;
            E0.block (i1, 0, i2 - i1 + 1,
                      E0.cols ())
                = fctr.transpose ()
                  * E0.block (i1, 0, i2 - i1 + 1,
                              E0.cols ());

            E1.block (0, i1, E1.rows (),
                      i2 - i1 + 1) *= fctr;
            E1.block (i1, 0, i2 - i1 + 1,
                      E1.cols ())
                = fctr.transpose ()
                  * E1.block (i1, 0, i2 - i1 + 1,
                              E1.cols ());

            E2.block (0, i1, E2.rows (),
                      i2 - i1 + 1) *= fctr;
            E2.block (i1, 0, i2 - i1 + 1,
                      E2.cols ())
                = fctr.transpose ()
                  * E2.block (i1, 0, i2 - i1 + 1,
                              E2.cols ());
        }

    E0 = (E0.array ().abs () < 1e-12)
             .select (Eigen::MatrixXd::Zero (
                          E0.rows (),E0.cols ()),E0);

    return { E0, E1, E2 };
}


Eigen::MatrixXd
computeInverseNormE0M0 (std::vector<NormE0Data> &NormE0Array, Eigen::MatrixXd &v,
                  int nd)
{
    for (auto &ii2 : NormE0Array)
        {
            int i1 = ii2.idx (0);
            int i2 = ii2.idx (1);
            Eigen::MatrixXd fctr = ii2.fctr;

            v.block (i1, 0, i2 - i1 + 1,
                     v.cols ())
                = fctr * v.block (i1, 0, i2 - i1 + 1,
                                  v.cols ());

            v.block (nd + i1, 0, i2 - i1 + 1,
                     v.cols ())
                = fctr.transpose ().fullPivLu ().solve (
                    v.block (nd + i1, 0,
                             i2 - i1 + 1, v.cols ()));
        }

    return v;
}



double calculateConditionNumber(const Eigen::MatrixXd &matrix)
{
    return 1.0 / matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
        .singularValues()
        .minCoeff();
}

Eigen::MatrixXd
matrixZMethod (Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2,
         int ndn)
{
    int nd = E0.rows ();
    Eigen::MatrixXd z = Eigen::MatrixXd::Zero (2 * nd, 2 * nd);
    Eigen::MatrixXd E0inv = E0.inverse ();

    z.block (0, nd, nd, nd) = -E0inv;
    z.block (nd, nd, nd, nd) = -E1 * E0inv;
    z.block (nd, 0, nd, nd)
        = -E2 + (-z.block (nd, nd, nd, nd)
                     * E1.transpose () - 1e-12 * E0);
    z.block (0, 0, nd, nd) =
        -(z.block (nd, nd, nd, nd))
             .transpose ();

    // calculate condition number
    double cond = calculateConditionNumber(z);
    std::cout << "COND Nr. Z-matrix: " << cond << std::endl;

    return z;
}