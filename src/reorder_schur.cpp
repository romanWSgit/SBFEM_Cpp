//
// Created by Roman Wallner- Silberhuber on 13.07.23.
//

#include "reorder_schur.h"

/**
 * @class SchurData
 * @brief Class for storing and processing data obtained from a real Schur decomposition.
 *
 * The SchurData class provides a way to store and process data obtained from a real Schur decomposition of a matrix.
 * It stores the matrix U, which represents the Schur vectors, the matrix T, which represents the Schur form,
 * the diagonal elements of T, and the eigenvalues of T.
 */
SchurData::SchurData (const Eigen::RealSchur<Eigen::MatrixXd> &A)
    :
      U (A.matrixU ()),
      S (A.matrixT ()),
      diagonalElements (S.diagonal ()),
      maxDiagonalElement (diagonalElements.maxCoeff ()),
      minDiagonalElement (diagonalElements.minCoeff ()),
      eigenvaluesOfS (
          (Eigen::EigenSolver<Eigen::MatrixXd> (S)).eigenvalues ())
{
}

// SchurData::SchurData (const Eigen::RealSchur<Eigen::MatrixXd> &A)
//{
//     this->U = A.matrixU ();
//     this->S = A.matrixT ();
//
//     Eigen::EigenSolver<Eigen::MatrixXd> solver (S);
//     const Eigen::VectorXcd &eigenvalues = solver.eigenvalues ();
//     this->eigenvaluesOfS = eigenvalues;
//     this->diagonalElements = S.diagonal ();
//     this->maxDiagonalElement = diagonalElements.maxCoeff ();
//     this->minDiagonalElement = diagonalElements.minCoeff ();
// }

SchurData::SchurData (const Eigen::MatrixXd &Q, const Eigen::MatrixXd &S)
{
    this->U = Q;
    this->S = S;

    Eigen::EigenSolver<Eigen::MatrixXd> solver (S);
    const Eigen::VectorXcd &eigenvalues = solver.eigenvalues ();
    this->eigenvaluesOfS = eigenvalues;
    this->diagonalElements = S.diagonal ();
    this->maxDiagonalElement = diagonalElements.maxCoeff ();
    this->minDiagonalElement = diagonalElements.minCoeff ();
}

SchurData::~SchurData ()
{
    std::cout << "Schur data got destroyed" << std::endl;
}

void
SchurData::getMinMaxEv ()
{
    Eigen::VectorXcd eV = this->eigenvaluesOfS;
    double maxReal = eV[0].real ();
    int maxIndex = 0;
    for (int i = 1; i < eV.size (); i++)
        {
            if (eV[i].real () > maxReal)
                {
                    maxReal = eV[i].real ();
                    maxIndex = i;
                    this->maxEV = eV[i];
                }
        }

    double minReal = eV[0].real ();
    int minIndex = 0;
    for (int i = 1; i < eV.size (); i++)
        {
            if (eV[i].real () < minReal)
                {
                    minReal = eV[i].real ();
                    minIndex = i;
                    this->minEV = eV[i];
                }
        }
}

// std::vector<int> swaplist(Eigen::VectorXcd p, std::vector<int> s,
// std::complex<double> z, double b)
//{
//     int n = p.size();
//     int k = 0;
//     int srtd = 0;  // Number of sorted eigenvalues.
//     bool fini = false;
//     std::vector<int> v;
//
//     // Compute block sizes.
//     std::vector<int> q(s.size() - 1);
//     std::transform(s.begin()+1, s.end(), s.begin(), q.begin(),
//     std::minus<int>());
//
//     while (!fini) {
//         auto [_, j] = select(p.segment(k, n - k), z);// Using Eigen's
//         segment for slicing
//
//         auto temp_p = p[j + k];
//         auto temp_q = q[j + k];
//         std::cout << "n: " << n << " j: " << j <<  " k: " << k << std::endl;
//         // insert this block at position k and remove it from where it was
//         taken. p.segment(j + k, n - j - k) = p.segment(j + k + 1, n - j - k
//         - 1); p[k] = temp_p;
//
//         // Similar for the block-sizes
//         q.erase(q.begin() + j + k);
//         q.insert(q.begin() + k, temp_q);
//
//         // Update the list of swaps for this block
//         for (int i = j + k - 1; i >= k; i--)
//         {
//             v.push_back(i);
//         }
//
//         // Update the number of sorted eigenvalues
//         srtd += q[k];
//         k++;
//
//         fini = (k >= n - 1 || k == -b || srtd == b || (srtd == b + 1 && b !=
//         0));
//     }
//     return v;
// }

std::vector<int>
swaplist (Eigen::VectorXcd p, std::vector<int> s, std::complex<double> z,
          double b)
{
    Eigen::VectorXcd p_orig = p;
    int n = p.size ();
    int k = 0;
    std::vector<int> v;
    int srtd = 0;
    std::vector<int> q (s.size () - 1);
    std::adjacent_difference (s.begin (), s.end (), q.begin ());
    q.erase (q.begin ()); // The first value is incorrect after
                          // std::adjacent_difference
    bool fini = false;
    while (!fini)
        {
            auto [_, j] = select (p.segment (k, n - k),
                                  z); // Using Eigen's segment for slicing
            std::complex<double> p_j = p[k + j];
            p.segment (k + 1, n - k - 1) = p.segment (k, n - k - 1);
            p[k] = p_j;
            p.conservativeResize (n - 1);
            int q_j = q[k + j];
            q.insert (q.begin () + k, q_j);
            q.erase (q.begin () + k + j + 1);
            for (int i = j + k; i >= k; --i)
                {
                    v.push_back (i);
                }
            srtd += q[k];
            k += 1;
            fini = (k >= n - 1) || (k == -b) || (srtd == b)
                   || (srtd == b + 1 && b != 0);
        }
    return v;
}

std::pair<double, Eigen::Index>
select (const Eigen::VectorXcd &p, std::complex<double> z)
{
    if (std::isinf (std::abs (z)))
        {
            Eigen::Index pos;
            double max_val = p.cwiseAbs ().maxCoeff (&pos);
            return std::make_pair (-max_val, pos);
        }
    else
        {
            std::complex<double> y
                = std::complex<double> (z.real (), std::abs (z.imag ()));
            Eigen::VectorXd delta = (p.array () - y).cwiseAbs ();
            Eigen::Index pos;
            double min_val = delta.minCoeff (&pos);
            return std::make_pair (min_val, pos);
        }
}

Eigen::Matrix2d
rot (const Eigen::Matrix2d &X)
{
    double c = 1.0;
    double s = 0.0;

    if (X (0, 0) != X (1, 1))
        {
            double tau = (X (0, 1) + X (1, 0)) / (X (0, 0) - X (1, 1));
            double off = std::sqrt (std::pow (tau, 2) + 1);
            std::vector<double> v = { tau - off, tau + off };
            int w = (std::abs (v[0]) < std::abs (v[1])) ? 0 : 1;
            c = 1.0 / std::sqrt (1.0 + std::pow (v[w], 2));
            s = v[w] * c;
        }

    Eigen::Matrix2d Q;
    Q << c, -s, s, c;
    return Q;
}

/**
 * Applies a Givens rotation such that the two-by-two diagonal block of S
 * situated at diagonal positions v[0], v[1] is in standardized form.
 *
 * @param schurUS Instance of Schur_U_S structure (copied by ref...).
 * @param v Index vector.
 *
 * @return Schur_U_S Modified Schur_U_S structure.
 */
SchurData
SchurData::normalizeInPlace (SchurData &schurUS, const Eigen::VectorXi &v)
{
    // Assuming rot function is defined elsewhere and returns 2x2 rotation
    // matrix
    Eigen::MatrixXd Q = rot (schurUS.S.block<2, 2> (v[0], v[1]));

    // Apply Givens rotation
    schurUS.S (Eigen::all, v) = schurUS.S (Eigen::all, v) * Q;
    schurUS.S (v, Eigen::all) = Q.transpose () * schurUS.S (v, Eigen::all);
    schurUS.U (Eigen::all, v) = schurUS.U (Eigen::all, v) * Q;

    return schurUS;
}

/**
 * Applies a Givens rotation such that the two-by-two diagonal block of S
 * situated at diagonal positions v[0], v[1] is in standardized form.
 *
 * @param schurUS Instance of Schur_U_S structure.
 * @param v Index vector.
 * @return Schur_U_S Modified Schur_U_S structure.
 */
SchurData
normalize (SchurData schurUS, Eigen::VectorXi v)
{
    // Assuming rot function is defined elsewhere and returns 2x2 rotation
    // matrix
    Eigen::MatrixXd Q = rot (schurUS.S.block<2, 2> (v[0], v[1]));

    // Apply Givens rotation
    schurUS.S (Eigen::all, v) = schurUS.S (Eigen::all, v) * Q;
    schurUS.S (v, Eigen::all) = Q.transpose () * schurUS.S (v, Eigen::all);
    schurUS.U (Eigen::all, v) = schurUS.U (Eigen::all, v) * Q;

    return schurUS;
}

int
select_pos (std::vector<std::complex<double> > p, double z)
{
    double y;
    std::vector<double> delta;
    int pos;

    if (std::isinf (std::abs (z)))
        {
            double max = -std::numeric_limits<double>::infinity ();
            pos = 0;
            for (int i = 0; i < p.size (); ++i)
                {
                    double abs_value
                        = std::sqrt (std::pow (p[i].real (), 2)
                                     + std::pow (p[i].imag (), 2));
                    if (abs_value > max)
                        {
                            max = abs_value;
                            pos = i;
                        }
                }
            return pos;
        }
    else
        {
            y = z;
            for (auto &val : p)
                {
                    delta.push_back (std::sqrt (std::pow (val.real () - y, 2)
                                                + std::pow (val.imag (), 2)));
                }

            double min = std::numeric_limits<double>::infinity ();
            pos = 0;
            for (int i = 0; i < delta.size (); ++i)
                {
                    if (delta[i] < min)
                        {
                            min = delta[i];
                            pos = i;
                        }
                }
            return pos;
        }
}

// std::vector<int> swaplist(std::vector<std::complex<double>> p,
// std::vector<int> s, std::complex<double> z, int b)
//{
//     int n = p.size();
//     std::vector<int> v;
//     int srtd = 0;  // Number of sorted eigenvalues.
//     std::vector<int> q(n);
//     std::partial_sum(s.begin(), s.end(), q.begin(), std::minus<int>()); //
//     Compute block sizes.
//
//     bool fini = false;
//     int k = 0;
//
//     while (!fini) {
//         int j = select_pos(p, z.real());  // Determine which block will go
//         to position k
//
//         // insert this block at position k,
//         std::rotate(p.begin() + k, p.begin() + k + j, p.end());
//
//         // and remove it from where it was taken.
//         p.erase(p.begin() + j + k + 1);
//
//         // Similar for the block-sizes
//         std::rotate(q.begin() + k, q.begin() + k + j, q.end());
//
//         // Update the list of swaps for this block
//         for (int i = j + k; i > k; --i) {
//             v.push_back(i);
//         }
//
//         // Update the number of sorted eigenvalues
//         srtd += q[k];
//         k++;
//         fini = (k >= n - 1 || k == -b || srtd == b || (srtd == b + 1 && b !=
//         0));
//     }
//     return v;
// }

void
swap (SchurData schurData, Eigen::VectorXi v, Eigen::VectorXi w)
{
    int p = v.size ();
    int q = w.size ();
    Eigen::MatrixXd Ip = Eigen::MatrixXd::Identity (p, p);
    Eigen::MatrixXd Iq = Eigen::MatrixXd::Identity (q, q);
    Eigen::VectorXd r (q);
    for (int j = 0; j < q; ++j)
        {
            r (j) = schurData.S (v[j], w[j]);
        }
    Eigen::MatrixXd K = Eigen::kroneckerProduct (
                            Iq, schurData.S (v, Eigen::all)
                                    .block (0, v (0), v.size (), v.size ()))
                        - Eigen::kroneckerProduct (
                            schurData.S (w, Eigen::all)
                                .block (0, w (0), w.size (), w.size ())
                                .transpose (),
                            Ip);
    Eigen::FullPivLU<Eigen::MatrixXd> lu (
        K); // LU-decomposition of this matrix.
    double e = lu.maxPivot ();
    Eigen::VectorXd sigp = Eigen::VectorXd::LinSpaced (p * q, 0, p * q - 1);
    for (int k = 0; k < p * q - 1;
         ++k) // Implement permutation P of the LU-decomposition PAQ=LU ...
        std::swap (sigp (k), sigp (lu.permutationP ().indices () (k)));
    r = e * r (sigp.array ().cast<int> ());
    Eigen::VectorXd x = lu.solve (r); // and solve the two triangular systems.
    Eigen::VectorXd sigq = Eigen::VectorXd::LinSpaced (p * q, 0, p * q - 1);
    for (int k = 0; k < p * q - 1;
         ++k) // Implement permutation Q of the LU-decomposition PAQ=LU ...
        std::swap (sigq (k), sigq (lu.permutationQ ().indices () (k)));
    x (sigq.array ().cast<int> ()) = x;
    Eigen::MatrixXd X (p, q);
    for (int j = 0; j < q; ++j) // De-vectorize the solution back to a block,
                                // or, quit Kronecker formulation.
        X.col (j) = x.segment (j * p, p);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr (
        (-X).rowwise ()
            .reverse ()
            .colwise ()
            .reverse ()); // Householder QR-decomposition of X.
    Eigen::MatrixXd Q = qr.householderQ ();
    Eigen::VectorXi vw = v;
    vw.conservativeResize (v.size () + w.size ());
    vw.tail (w.size ()) = w;

    schurData.S.block (0, vw (0), schurData.S.rows (), vw.size ())
        = schurData.S.block (0, vw (0), schurData.S.rows (), vw.size ())
          * Q; // Perform the actual swap by left- and right-multiplication of
               // S by Q,
    schurData.S.block (vw (0), 0, vw.size (), schurData.S.cols ())
        = Q.transpose ()
          * schurData.S.block (vw (0), 0, vw.size (), schurData.S.cols ());
    schurData.U.block (0, vw (0), schurData.U.rows (), vw.size ())
        = schurData.U.block (0, vw (0), schurData.U.rows (), vw.size ())
          * Q; // and, right-multiplication of U by Q
}

void
check_block_triangular (Eigen::MatrixXd &R, double eps)
{
    int rows = R.rows ();
    int cols = R.cols ();

    for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
                {
                    if (i < j + 2 && std::abs (R (i, j)) > 100 * eps)
                        {
                            throw std::invalid_argument (
                                "R is not block-triangular");
                        }
                }
        }
}

Eigen::VectorXd
SRSchur (SchurData schurData, double z, int b)
{

    Eigen::VectorXd ap;
    int cnt = 1;
    std::vector<int> result2;
    std::vector<int> result3;
    std::vector<int> r;
    std::vector<int> s (schurData.S.rows () + 1);
    iota (s.begin (), s.end (),
          1); // Fill with increasing series starting from 1
    std::vector<std::complex<double> > p (s.size () - 1);

    // Getting the subdiagonal elements

    // R is Eigen::MatrixXd
    Eigen::VectorXd subDiag = schurData.S.diagonal (-1).cwiseAbs ();
    for (int i = 0; i < subDiag.size (); i++)
        {
            if (subDiag[i] > 100 * std::numeric_limits<double>::epsilon ())
                {
                    r.push_back (i);
                }
        }

    //    for (int i = 1; i < schurData.S.rows(); ++i)
    //    {
    //        if (abs(schurData.S(i, i - 1)) > 100 *
    //        std::numeric_limits<double>::epsilon())
    //            r.push_back(i);
    //    }

    // Removing the identified subdiagonal elements
    for (auto &ri : r)
        s.erase (s.begin () + ri);

    for (int k = 0; k < s.size () - 1; ++k)
        {
            int sk = s[k];
            if (s[k + 1] - sk == 2)
                {
                    int size = s[k + 1] - sk;
                    std::vector<int> result1 (size);
                    int start = s[k];
                    int end = s[k + 1] - 1;
                    for (int i = 0; i < size; ++i)
                        {
                            result1[i] = start + i;
                        }
                    Eigen::Map<Eigen::VectorXi> result1_eigen (
                        result1.data (), result1.size ());
                    normalize (schurData, result1_eigen);
                    std::complex<double> p1 = schurData.S (sk - 1, sk - 1);
                    std::complex<double> comp (schurData.S (sk, sk - 1)
                                                   * schurData.S (sk - 1, sk),
                                               0.0);
                    std::complex<double> p2 = std::sqrt (comp);
                    p[k] = p1 + p2;
                }
            else
                {
                    std::complex<double> p1 = schurData.S (s[k] - 1, s[k] - 1);
                    p[k] = p1;
                }
        }

    Eigen::Map<Eigen::VectorXcd> p_eigen (p.data (), p.size ());
    std::vector<int> swapEXP
        = swaplist (p_eigen, s, std::complex<double> (z, 0.0), b);
    ap = Eigen::VectorXd::Zero (swapEXP.size ());
    for (int k = 0; k < swapEXP.size (); k++)
        {
            std::vector<int> v (s.begin () + swapEXP[k],
                                s.begin () + swapEXP[k + 1]);
            std::vector<int> w (s.begin () + swapEXP[k + 1],
                                s.begin () + swapEXP[k + 2]);
            Eigen::VectorXi vEigen
                = Eigen::Map<Eigen::VectorXi> (v.data (), v.size ());
            Eigen::VectorXi wEigen
                = Eigen::Map<Eigen::VectorXi> (w.data (), w.size ());
            double nrA = (schurData.S
                              .block (v.front () - 1, v.front () - 1,
                                      v.size (), v.size ())
                              .array ()
                              .abs ())
                             .maxCoeff ();
            swap (schurData, vEigen, wEigen);
            s[swapEXP[k] + 1]
                = s[swapEXP[k]] + s[swapEXP[k] + 2] - s[swapEXP[k] + 1];
            v = std::vector<int> (s.begin () + swapEXP[k],
                                  s.begin () + swapEXP[k + 1]);
            w = std::vector<int> (s.begin () + swapEXP[k + 1],
                                  s.begin () + swapEXP[k + 2]);
            if (v.size () == 2)
                {
                    result2.push_back (v.front ());
                    result2.push_back (v.back ());
                    Eigen::Map<Eigen::VectorXi> result2_eigen (
                        result2.data (), result2.size ());
                    normalize (schurData, result2_eigen);
                }
            if (w.size () == 2)
                {
                    result3.push_back (w.front ());
                    result3.push_back (w.back ());
                    Eigen::Map<Eigen::VectorXi> result3_eigen (
                        result3.data (), result3.size ());
                    normalize (schurData, result3_eigen);
                }
            ap[k] = (schurData.S
                         .block (w.front () - 1, v.front () - 1, w.size (),
                                 v.size ())
                         .array ()
                         .abs ())
                        .maxCoeff ()
                    / (10 * std::numeric_limits<double>::epsilon () * nrA);
        }

    schurData.S -= schurData.S.triangularView<Eigen::StrictlyLower> ();
    for (int j = 2; j < s.size () - 1; ++j)
        schurData.S (s[j], s[j] - 1) = 0;
    return ap;
}
