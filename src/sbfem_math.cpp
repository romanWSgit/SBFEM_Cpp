//
// Created by Roman Wallner- Silberhuber on 20.05.23.
//

#include "sbfem_math.h"
#include "ShapeFunctionType.h"

Eigen::VectorXd evenDistributedLagrangianPoints(int polyOrd)
{
    Eigen::VectorXd points(polyOrd + 1);
    for (int i = 0; i <= polyOrd; ++i)
    {
        points(i) = -1 + (i * 2.0) / polyOrd;
    }
    return points;
}

double lagrange(double x, int i, Eigen::VectorXd xm)
{
    int n = xm.size() - 1;
    double y = 1.0;
    for (int j = 0; j <= n; j++)
    {
        if (i != j)
            y *= (x - xm(j)) / (xm(i) - xm(j));
    }
    return y;
}

double lagrange_diff(double x, int j, Eigen::VectorXd xm)
{
    int n = xm.size() - 1;
    double y = 0.0;
    for (int l = 0; l <= n; l++)
    {
        if (l != j)
        {
            double k = 1.0 / (xm(j) - xm(l));
            for (int m = 0; m <= n; m++)
            {
                if (m != j && m != l)
                    k *= (x - xm(m)) / (xm(j) - xm(m));
            }
            y += k;
        }
    }
    return y;
}

double legendre_poly_integrated(double x, int p)
{
    switch (p)
    {
    case -1:
        return 0.5 - 0.5 * x;
    case 0:
        return 0.5 + 0.5 * x;
    case 1:
        return -1.0 + std::pow(x, 2);
    case 2:
        return 2.0 * x * (-1.0 + std::pow(x, 2));
    case 3:
        return 0.75 - 4.5 * std::pow(x, 2) + 3.75 * std::pow(x, 4);
    case 4:
        return 3.0 * x - 10.0 * std::pow(x, 3) + 7.0 * std::pow(x, 5);
    case 5:
        return -0.625 + 9.375 * std::pow(x, 2) - 21.875 * std::pow(x, 4) +
               13.125 * std::pow(x, 6);
    case 6:
        return -3.75 * x + 26.25 * std::pow(x, 3) - 47.25 * std::pow(x, 5) +
               24.75 * std::pow(x, 7);
    case 7:
        return 0.546875 - 15.3125 * std::pow(x, 2) + 68.90625 * std::pow(x, 4) -
               101.0625 * std::pow(x, 6) + 46.921875 * std::pow(x, 8);
    case 8:
        return 4.375 * x - 52.5 * std::pow(x, 3) + 173.25 * std::pow(x, 5) -
               214.5 * std::pow(x, 7) + 89.375 * std::pow(x, 9);
    case 9:
        return -0.4921875 + 22.1484375 * std::pow(x, 2) -
               162.421875 * std::pow(x, 4) + 422.296875 * std::pow(x, 6) -
               452.4609375 * std::pow(x, 8) + 170.9296875 * std::pow(x, 10);
    case 10:
        return -4.921875 * x + 90.234375 * std::pow(x, 3) -
               469.21875 * std::pow(x, 5) + 1005.46875 * std::pow(x, 7) -
               949.609375 * std::pow(x, 9) + 328.046875 * std::pow(x, 11);
    case 11:
        return 0.451171875 - 29.77734375 * std::pow(x, 2) +
               322.587890625 * std::pow(x, 4) - 1290.3515625 * std::pow(x, 6) +
               2350.283203125 * std::pow(x, 8) -
               1984.68359375 * std::pow(x, 10) +
               631.490234375 * std::pow(x, 12);
    case 12:
        return 5.4140625 * x - 140.765625 * std::pow(x, 3) +
               1055.7421875 * std::pow(x, 5) - 3418.59375 * std::pow(x, 7) +
               5412.7734375 * std::pow(x, 9) - 4133.390625 * std::pow(x, 11) +
               1218.8203125 * std::pow(x, 13);
    case 13:
        return -0.4189453125 + 38.1240234375 * std::pow(x, 2) -
               571.8603515625 * std::pow(x, 4) +
               3240.5419921875 * std::pow(x, 6) -
               8795.7568359375 * std::pow(x, 8) +
               12314.0595703125 * std::pow(x, 10) -
               8582.5263671875 * std::pow(x, 12) +
               2357.8369140625 * std::pow(x, 14);
    case 14:
        return -5.865234375 * x + 205.283203125 * std::pow(x, 3) -
               2093.888671875 * std::pow(x, 5) +
               9472.353515625 * std::pow(x, 7) -
               22102.158203125 * std::pow(x, 9) +
               27728.162109375 * std::pow(x, 11) -
               17774.462890625 * std::pow(x, 13) +
               4570.576171875 * std::pow(x, 15);
    case 15:
        return 0.39276123046875 - 47.13134765625 * std::pow(x, 2) +
               934.771728515625 * std::pow(x, 4) -
               7104.26513671875 * std::pow(x, 6) +
               26640.994262695312 * std::pow(x, 8) -
               54466.03271484375 * std::pow(x, 10) +
               61893.218994140625 * std::pow(x, 12) -
               36727.84423828125 * std::pow(x, 14) +
               8875.895690917969 * std::pow(x, 16);
    case 16:
        return 6.2841796875 * x - 284.8828125 * std::pow(x, 3) +
               3788.94140625 * std::pow(x, 5) - 22733.6484375 * std::pow(x, 7) +
               72621.376953125 * std::pow(x, 9) -
               132038.8671875 * std::pow(x, 11) +
               137117.28515625 * std::pow(x, 13) -
               75740.9765625 * std::pow(x, 15) +
               17264.4873046875 * std::pow(x, 17);
    case 17:
        return 0.370941162109375 - 56.753997802734375 * std::pow(x, 2) +
               1437.7679443359375 * std::pow(x, 4) -
               14090.125854492188 * std::pow(x, 6) +
               69444.19171142578 * std::pow(x, 8) -
               192900.53253173828 * std::pow(x, 10) +
               315655.4168701172 * std::pow(x, 12) -
               301780.45349121094 * std::pow(x, 14) +
               155919.90097045898 * std::pow(x, 16) -
               33629.78256225586 * std::pow(x, 18);
    case 18:
        return -6.67694091796875 * x + 380.58563232421875 * std::pow(x, 3) -
               6393.838623046875 * std::pow(x, 5) +
               49019.429443359375 * std::pow(x, 7) -
               204247.62268066406 * std::pow(x, 9) +
               501335.07385253906 * std::pow(x, 11) -
               745575.2380371094 * std::pow(x, 13) +
               660366.6394042969 * std::pow(x, 15) -
               320472.0455932617 * std::pow(x, 17) +
               65593.69354248047 * std::pow(x, 19);
    case 19:
        return 0.35239410400390625 - 66.95487976074219 * std::pow(x, 2) +
               2109.078712463379 * std::pow(x, 4) -
               25871.36553955078 * std::pow(x, 6) +
               161696.03462219238 * std::pow(x, 8) -
               582105.7246398926 * std::pow(x, 10) +
               1.2788686374664307e6 * std::pow(x, 12) -
               1.7426341873168945e6 * std::pow(x, 14) +
               1.437673204536438e6 * std::pow(x, 16) -
               657758.9824676514 * std::pow(x, 18) +
               128089.90711212158 * std::pow(x, 20);
    case 20:
        return 7.047882080078125 * x - 493.35174560546875 * std::pow(x, 3) +
               10212.381134033203 * std::pow(x, 5) -
               97260.77270507812 * std::pow(x, 7) +
               510619.05670166016 * std::pow(x, 9) -
               1.6154130157470703e6 * std::pow(x, 11) +
               3.210115608215332e6 * std::pow(x, 13) -
               4.035573907470703e6 * std::pow(x, 15) +
               3.1157004432678223e6 * std::pow(x, 17) -
               1.348314811706543e6 * std::pow(x, 19) +
               250401.32217407227 * std::pow(x, 21);
    default:
        std::cerr << "Case not implemented.\n";
        return 0.0;
    }
}

double legendre_poly_integrated_diff(double x, int p)
{
    switch (p)
    {
    case -1:
        return -0.5;
    case 0:
        return 0.5;
    case 1:
        return 2. * x;
    case 2:
        return -2. + 6. * std::pow(x, 2);
    case 3:
        return -9. * x + 15. * std::pow(x, 3);
    case 4:
        return 3. - 30. * std::pow(x, 2) + 35. * std::pow(x, 4);
    case 5:
        return 18.75 * x - 87.5 * std::pow(x, 3) + 78.75 * std::pow(x, 5);
    case 6:
        return -3.75 + 78.75 * std::pow(x, 2) - 236.25 * std::pow(x, 4) +
               173.25 * std::pow(x, 6);
    case 7:
        return -30.625 * x + 275.625 * std::pow(x, 3) -
               606.375 * std::pow(x, 5) + 375.375 * std::pow(x, 7);
    case 8:
        return 4.375 - 157.5 * std::pow(x, 2) + 866.25 * std::pow(x, 4) -
               1501.5 * std::pow(x, 6) + 804.375 * std::pow(x, 8);
    case 9:
        return 44.296875 * x - 649.6875 * std::pow(x, 3) +
               2533.78125 * std::pow(x, 5) - 3619.6875 * std::pow(x, 7) +
               1709.296875 * std::pow(x, 9);
    case 10:
        return -4.921875 + 270.703125 * std::pow(x, 2) -
               2346.09375 * std::pow(x, 4) + 7038.28125 * std::pow(x, 6) -
               8546.484375 * std::pow(x, 8) + 3608.515625 * std::pow(x, 10);
    case 11:
        return -59.5546875 * x + 1290.3515625 * std::pow(x, 3) -
               7742.109375 * std::pow(x, 5) + 18802.265625 * std::pow(x, 7) -
               19846.8359375 * std::pow(x, 9) + 7577.8828125 * std::pow(x, 11);
    case 12:
        return 5.4140625 - 422.296875 * std::pow(x, 2) +
               5278.7109375 * std::pow(x, 4) - 23930.15625 * std::pow(x, 6) +
               48714.9609375 * std::pow(x, 8) - 45467.296875 * std::pow(x, 10) +
               15844.6640625 * std::pow(x, 12);
    case 13:
        return 76.248046875 * x - 2287.44140625 * std::pow(x, 3) +
               19443.251953125 * std::pow(x, 5) -
               70366.0546875 * std::pow(x, 7) +
               123140.595703125 * std::pow(x, 9) -
               102990.31640625 * std::pow(x, 11) +
               33009.716796875 * std::pow(x, 13);
    case 14:
        return -5.865234375 + 615.849609375 * std::pow(x, 2) -
               10469.443359375 * std::pow(x, 4) +
               66306.474609375 * std::pow(x, 6) -
               198919.423828125 * std::pow(x, 8) +
               305009.783203125 * std::pow(x, 10) -
               231068.017578125 * std::pow(x, 12) +
               68558.642578125 * std::pow(x, 14);
    case 15:
        return -94.2626953125 * x + 3739.0869140625 * std::pow(x, 3) -
               42625.5908203125 * std::pow(x, 5) +
               213127.9541015625 * std::pow(x, 7) -
               544660.3271484375 * std::pow(x, 9) +
               742718.6279296875 * std::pow(x, 11) -
               514189.8193359375 * std::pow(x, 13) +
               142014.3310546875 * std::pow(x, 15);
    case 16:
        return 6.2841796875 - 854.6484375 * std::pow(x, 2) +
               18944.70703125 * std::pow(x, 4) -
               159135.5390625 * std::pow(x, 6) +
               653592.392578125 * std::pow(x, 8) -
               1.4524275390625e6 * std::pow(x, 10) +
               1.78252470703125e6 * std::pow(x, 12) -
               1.1361146484375e6 * std::pow(x, 14) +
               293496.2841796875 * std::pow(x, 16);
    case 17:
        return 113.50799560546875 * x - 5751.07177734375 * std::pow(x, 3) +
               84540.75512695312 * std::pow(x, 5) -
               555553.5336914062 * std::pow(x, 7) +
               1.9290053253173828e6 * std::pow(x, 9) -
               3.7878650024414068e6 * std::pow(x, 11) +
               4.224926348876953e6 * std::pow(x, 13) -
               2.4947184155273438e6 * std::pow(x, 15) +
               605336.0861206055 * std::pow(x, 17);
    case 18:
        return -6.67694091796875 + 1141.7568969726562 * std::pow(x, 2) -
               31969.193115234375 * std::pow(x, 4) +
               343136.0061035156 * std::pow(x, 6) -
               1.8382286041259766e6 * std::pow(x, 8) +
               5.51468581237793e6 * std::pow(x, 10) -
               9.692478094482422e6 * std::pow(x, 12) +
               9.905499591064453e6 * std::pow(x, 14) -
               5.448024775085449e6 * std::pow(x, 16) +
               1.246280177307129e6 * std::pow(x, 18);
    case 19:
        return 133.90975952148438 * x - 8436.314849853516 * std::pow(x, 3) +
               155228.1932373047 * std::pow(x, 5) -
               1.293568276977539e6 * std::pow(x, 7) +
               5.821057246398926e6 * std::pow(x, 9) -
               1.534642364959717e7 * std::pow(x, 11) +
               2.4396878622436523e7 * std::pow(x, 13) -
               2.300277127258301e7 * std::pow(x, 15) +
               1.1839661684417725e7 * std::pow(x, 17) -
               2.5617981422424316e6 * std::pow(x, 19);
    case 20:
        return -7.047882080078125 + 1480.0552368164062 * std::pow(x, 2) -
               51061.905670166016 * std::pow(x, 4) +
               680825.4089355469 * std::pow(x, 6) -
               4.5955715103149414e6 * std::pow(x, 8) +
               1.7769543173217774e7 * std::pow(x, 10) -
               4.1731502906799316e7 * std::pow(x, 12) +
               6.053360861206055e7 * std::pow(x, 14) -
               5.296690753555298e7 * std::pow(x, 16) +
               2.5617981422424313e7 * std::pow(x, 18) -
               5.258427765655518e6 * std::pow(x, 20);
    default:
        std::cerr << "Case not implemented.\n";
        return 0.0;
    }
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXd> shape_N(ShapeFunctionType shapeFct,
                                                     double eta, int poly_ord)
{
    if (shapeFct == ShapeFunctionType::STANDARD)
    {
        std::vector<double> n_vec;
        Eigen::MatrixXd n_mat(2, 0);
        for (int i = 0; i <= poly_ord; ++i)
        {
            double result =
                lagrange(eta, i, evenDistributedLagrangianPoints(poly_ord));
            n_vec.push_back(result);
            n_mat.conservativeResize(2, n_mat.cols() + 2);
            n_mat.col(n_mat.cols() - 2) << result, 0;
            n_mat.col(n_mat.cols() - 1) << 0, result;
        }
        return {Eigen::VectorXd::Map(n_vec.data(), n_vec.size()), n_mat};
    }
    else if (shapeFct == ShapeFunctionType::HIERARCHICAL)
    {
        std::vector<double> n_vec;
        Eigen::MatrixXd n_mat(2, 0);
        for (int i = -1; i <= poly_ord - 1; ++i)
        {
            double result = legendre_poly_integrated(eta,i);
            n_vec.push_back(result);
            n_mat.conservativeResize(2, n_mat.cols() + 2);
            n_mat.col(n_mat.cols() - 2) << result, 0;
            n_mat.col(n_mat.cols() - 1) << 0, result;
        }
        return {Eigen::VectorXd::Map(n_vec.data(), n_vec.size()), n_mat};
    }
    else
    {
        throw std::runtime_error("Invalid shape function type.");
    }
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXd> shape_dN(double eta, int poly_ord)
{
    std::vector<double> n_vec;
    Eigen::MatrixXd n_mat(2, 0);
    for (int i = 0; i <= poly_ord; i++)
    {
        double result =
            lagrange_diff(eta, i, evenDistributedLagrangianPoints(poly_ord));
        n_vec.push_back(result);
        n_mat.conservativeResize(2, n_mat.cols() + 2);
        n_mat.col(n_mat.cols() - 2) << result, 0;
        n_mat.col(n_mat.cols() - 1) << 0, result;
    }
    return {Eigen::VectorXd::Map(n_vec.data(), n_vec.size()), n_mat};
}
