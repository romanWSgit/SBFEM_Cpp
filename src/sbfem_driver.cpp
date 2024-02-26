//
// Created by Roman Wallner- Silberhuber on 20.02.24.
//
#import "sbfem_driver.h"

// Example function that takes a scalar and returns a 2D vector
Eigen::VectorXd exampleFunction(double x)
{
    Eigen::VectorXd result(2);
    result << sin(x), cos(x); // Just as an example
    return result;
}

// Wrapper to integrate a specific component of the vector returned by
// exampleFunction
double integrateComponent(double a, double b, int componentIndex)
{
    auto f = [componentIndex](double x) {
        Eigen::VectorXd v = exampleFunction(x);
        return v(componentIndex);
    };

    // Use Gauss-Kronrod quadrature for integration
    auto q = boost::math::quadrature::gauss_kronrod<double, 15>();
    // Instantiate Gauss-Kronrod with 15 points
    double result = q.integrate(f, a, b);
    return result;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, double, double> matrix_integrator(
    double a, double b, const std::function<Eigen::MatrixXd(double)> &matFunc,
    int rows, int cols)
{
    Eigen::MatrixXd integralResults(rows, cols);
    Eigen::MatrixXd errorEstimates(rows, cols);

    double maxError = 0;
    double minError = std::numeric_limits<double>::infinity();

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            auto elementFunc = [&, i, j](double x) { return matFunc(x)(i, j); };

            double error_estimate;
            double result =
                boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
                    elementFunc, a, b, 5, 1e-10, &error_estimate);

            integralResults(i, j) = result;
            errorEstimates(i, j) = error_estimate;

            if (error_estimate > maxError)
                maxError = error_estimate;
            if (error_estimate < minError)
                minError = error_estimate;
        }
    }

    return {integralResults, errorEstimates, maxError, minError};
}

double f_new(double t)
{
    return std::exp(-t * t / 2);
};

double integrateComponentKronrod(double a, double b, int componentIndex)
{
    auto f = [componentIndex](double x) {
        Eigen::VectorXd v = exampleFunction(x);
        return v(componentIndex);
    };
    double error_estimate;
    // Use Gauss-Kronrod quadrature for integration
    auto q = boost::math::quadrature::gauss_kronrod<double, 15>();
    // Instantiate Gauss-Kronrod with 15 points
    double result = q.integrate(f, a, b, 5, 1e-10, &error_estimate);
    std::cout << "Component " << componentIndex
              << " integration error estimate: " << error_estimate << std::endl;
    return result;
}

[[maybe_unused]] const int superElements = 1;

//void main_loop(const SuperElementJson &sE)
//{
//    for (int se = 1; se <= superElements; ++se)
//    {
//        int nodesPerElement = sE.getMSENodedim() eDim = nodesPerElement * ndn;
//        int nElements = sE.getMSEIelem();
//        for (int el = 1; el <= nElements; ++el)
//        {
//            Eigen::MatrixXd matrix(rows, cols); // MatrixXd for double precision
//            matrix.setZero();
//        }
//    }
//}

// For[i = 1, i <= superElements, i++,
// for (int se = 0; se <= superElements; ++se)
//{
//     for (int el = 0; el < ielem; ++ielem)
//     {
//         e0 = e0 + intF*Transpose[B1elIp] . DmatParams . B1elIp;
//         e1 = e1 + intF*Transpose[B2elIp] . DmatParams . B1elIp;
//         e2 = e2 + intF*Transpose[B2elIp] . DmatParams . B2elIp;
//         m0 = m0 + intF*Transpose[NmatIP] . NmatIP*\[Rho]Params;
//         f0 = f0 +
//              intForce*
//                  Transpose[NmatIP] . If[fCond[i, el] == True, flastIp, {0.,
//                  0.}];
//     }
//     // code here
// }

// ndn = If[piezoFlag == True, 3, 2]; (*number of dofs per node*)
//     nodesPerElement = polyOrd + 1;  (*number of nodes per element*)
//     addInt = 4; (*add to integration Points*)
//     nIntegrationPoints = polyOrd + 1 + addInt;
// nIntegrationPoints = ordGauss + addInt;
// intPointPos =
//     If[GLL, RootsLobatto[nIntegrationPoints],
//        GaussianQuadratureWeights[nIntegrationPoints, -1, 1][[All, 1]]];
// intPointWeights = If[GLL, Table[WeightsLobatto[
//                                     nIntegrationPoints, i], {i,
//                                     nIntegrationPoints}],
//                      GaussianQuadratureWeights[nIntegrationPoints, -1,
//                      1][[All, 2]]];
// CoordMat = {First /@ Partition[CoordVec, 2],
//             Last /@ Partition[CoordVec, 2]};
// NPlusDNMatrix[\[Eta]_] := {Nvec[\[Eta]], DNvec[\[Eta]]};
// eDim = nodesPerElement*ndn;
// aZeroDot =
//  Table[{Global`ax[i // N] -> 0, Global`ay[i // N] -> 0}, {i, 3,
//     polyOrd + 1}] // Flatten; E0Tx =
//  E1Tx = E2Tx =
//    M0Tx = F0Tx = E0Tx2 = E1Tx2 = E2Tx2 = M0Tx2 = F0Tx2 = {};
// For[i = 1, i <= superElements, i++,
//     E0TxElList =
//         E1TxElList =
//             E2TxElList =
//                 M0TxElList =
//                     F0TxElList =
//                         E0TxElList2 =
//                             E1TxElList2 = E2TxElList2 = M0TxElList2 =
//                             F0TxElList2 = {};
//     For[el = 1, el <= sE[[i, "ielem"]], el++,
//         coordMatEL = CoordMat /. sE[[i, "elementsC"]][[el]] // N;
//                                B1el[\[Eta]_] = B1[\[Eta]] /. sE[[i,
//                                "elementsC"]][[el]] // N;
//                                                             B2el[\[Eta]_] =
//                                                             B2[\[Eta]] /.
//                                                             sE[[i,
//                                                             "elementsC"]][[el]]
//                                                             // N;
//                                                                                          b1el[\[Eta]_] = b1[\[Eta]] /. sE[[i, "elementsC"]][[el]] // N;
//                                                                                                                       b2el[\[Eta]_] = b2[\[Eta]] /. sE[[i, "elementsC"]][[el]] // N;
//                                                                                                                                                    (*flastel[\[Eta]_]=flastP[i,el,\[Eta]]/.sE[[i,"elementsC"]][[el]]//
//                                                                                                                                                                                                N;*)
//                                                                                                                                                        flastel[\[Eta]_] =
//                                                                                                                                            fLast[i, el, \[Eta]] /. sE[[i, "elementsC"]][[el]] // N;
//                                                                                                                                                                  parEl = paramsL[el];
//         refEl = ref;
//         dCF = If[ndFlag == True, 1/(Er /. ref), 1.] /. parEl // N;
//                                                           DmatParams = Dmat
//                                                           /. parEl;
//    \[Rho]Params = \[Rho] /. parEl;
//         e0 = ConstantArray[0, {eDim, eDim}];
//         e1 = ConstantArray[0, {eDim, eDim}];
//         e2 = ConstantArray[0, {eDim, eDim}];
//         m0 = ConstantArray[0, {eDim, eDim}];
//         f0 = ConstantArray[0, {eDim}];
//         If[GLL,
//            e0n2 = ConstantArray[0, {eDim, eDim}];
//            e1n2 = ConstantArray[0, {eDim, eDim}];
//            e2n2 = ConstantArray[0, {eDim, eDim}];
//            m0n2 = ConstantArray[0, {eDim, eDim}];
//            f0n2 = ConstantArray[0, {eDim}];
//];
//         For[iP = 1, iP <= nIntegrationPoints, iP++,
//             If[GLL, kidx = Range[ndn * (iP - 1) + 1, ndn*iP];];
//     \[Eta]IP = intPointPos[[iP]];
//             jMatIp =
//                 NPlusDNMatrix[\[Eta]IP] . Transpose[coordMatEL] /. aZeroDot;
//             detJMatIP = Det[jMatIp];
//             invJMatIP = Inverse[jMatIp];
//             intF = intPointWeights[[iP]]*dCF*detJMatIP;
//             intForce =
//                 intPointWeights[[iP]]*dCF*
//                 Sqrt[g\[Xi][\[Eta]IP] . g\[Xi][\[Eta]IP]] /.
//                                                            sE[[i]]["elementsC"][[el]];
//             NmatIP =
//                 If[piezoFlag == True, NmatPiezo[\[Eta]IP], Nmat[\[Eta]IP]];
//             B1elIp = B1el[\[Eta]IP];
//             B2elIp = B2el[\[Eta]IP];
//             b1elIp = b1el[\[Eta]IP];
//             b2elIp = b2el[\[Eta]IP];
//             flastIp = flastel[\[Eta]IP];
//             mtmp = DmatParams . b1elIp;
//             If[GLL,
//                e0n2[[kidx, kidx]] = intF*Transpose[b1elIp] . mtmp;
//                e1n2[[All, kidx]] = intF*Transpose[B2elIp] . mtmp;
//                m0n2[[kidx, kidx]] = intF*\[Rho]Params*IdentityMatrix[ndn];
//                f0n2[[kidx]] =
//                    intForce*If[fCond[i, el] == True, flastIp, {0., 0.}];
//                e2n2 = e2n2 + intF*Transpose[B2elIp] . DmatParams . B2elIp;
//];
//
//             e0 = e0 + intF*Transpose[B1elIp] . DmatParams . B1elIp;
//             e1 = e1 + intF*Transpose[B2elIp] . DmatParams . B1elIp;
//             e2 = e2 + intF*Transpose[B2elIp] . DmatParams . B2elIp;
//             m0 = m0 + intF*Transpose[NmatIP] . NmatIP*\[Rho]Params;
//             f0 =
//                 f0 + intForce*
//                          Transpose[NmatIP] . If[fCond[i, el] == True,
//                          flastIp, {0., 0.}];
//];
//         AppendTo[E0TxElList, e0];
//         AppendTo[E1TxElList, e1];
//         AppendTo[E2TxElList, e2];
//         AppendTo[M0TxElList, m0];
//         AppendTo[F0TxElList, f0];
//         If[GLL,
//            AppendTo[E0TxElList2, e0n2];
//            AppendTo[E1TxElList2, e1n2];
//            AppendTo[E2TxElList2, e2n2];
//            AppendTo[M0TxElList2, m0n2];
//            AppendTo[F0TxElList2, f0n2];
//];
//];
//     AppendTo[E0Tx, E0TxElList];
//     AppendTo[E1Tx, E1TxElList];
//     AppendTo[E2Tx, E2TxElList];
//     AppendTo[M0Tx, M0TxElList];
//     AppendTo[F0Tx, F0TxElList];
//     If[GLL,
//        AppendTo[E0Tx2, E0TxElList2];
//        AppendTo[E1Tx2, E1TxElList2];
//        AppendTo[E2Tx2, E2TxElList2];
//        AppendTo[M0Tx2, M0TxElList2];
//        AppendTo[F0Tx2, F0TxElList2];
//];
//];