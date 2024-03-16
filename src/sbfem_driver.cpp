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

// just a test function
Eigen::MatrixXd myMatrixFunction(double x, int i)
{
    Eigen::MatrixXd result(2, 2); // Example: Create a 2x2 matrix
    result(0, 0) = cos(x) * i * x;
    result(0, 1) = x * x * i;
    result(1, 0) = pow(x, i);
    result(1, 1) = sin(x * i);
    return result;
}


void main_loop(const SuperElementJson &sE,
               std::vector<StructData> &material_data_list)
{
    omp_set_num_threads(10);     // Requests OpenMP to use 10 threads
    const int superElements = 1; // set here manually
    const int numberOfDegreesPerNode =
        2; // Note that this is limited to 2 at the moment
    const int nodesPerElement =
        sE.getMSEPolyOrd() + 1; // number of nodes per element ??? hierarchical
    // TODO: Thats a really good question! How does this work with
    // hierarchical???
    const int eDim = nodesPerElement * numberOfDegreesPerNode;
    const int rows = eDim;
    const int cols = eDim;

    // start clock
    auto start = std::chrono::high_resolution_clock::now();

    // here it is assumed that all all elements within a superelement
    // have the same poly_Ord and material data...etc...

    for (int se = 1; se <= superElements; ++se)
    {
        const Eigen::Vector2d centre = sE.getMSECentre();
        const int poly_Ord = sE.getMSEPolyOrd();
        const ShapeFunctionType shapeFunctionType = sE.getMSEShapeFct();

        std::string search_materialType = "Flat Disc";

        StructData material =
            getMaterialDataForType(search_materialType, material_data_list);

        auto material_data = findStructDataAndMaterial(
            material_data_list, "Flat Disc", "fictional");

        Eigen::MatrixXd Dmat = material_data->getDMatrix();
        std::cout << "Dmat: "
                  << "\n"
                  << Dmat << std::endl;
        const int elements = sE.getMSEIelem();

        std::vector<Eigen::MatrixXd> e0_list;
        std::vector<Eigen::MatrixXd> e1_list;
        std::vector<Eigen::MatrixXd> e2_list;

#pragma omp parallel
        {
            std::vector<Eigen::MatrixXd> e0_list_private;
            std::vector<Eigen::MatrixXd> e1_list_private;
            std::vector<Eigen::MatrixXd> e2_list_private;

#pragma omp for nowait // Each thread processes a subset of elements
            for (int el = 0; el < elements; ++el)
            {
                Eigen::VectorXd coord_vec = sE.getMSEElementsMC().row(el);

                Eigen::MatrixXd m0 = Eigen::MatrixXd::Zero(rows, cols);

                auto f_bind_e0 = std::bind(
                    ComputeEMatrixKernel, std::placeholders::_1, poly_Ord,
                    coord_vec, shapeFunctionType, centre, Dmat, B1);

                auto f_bind_e1 = std::bind(
                    ComputeMixedEMatrixKernel, std::placeholders::_1, poly_Ord,
                    coord_vec, shapeFunctionType, centre, Dmat);

                auto f_bind_e2 = std::bind(
                    ComputeEMatrixKernel, std::placeholders::_1, poly_Ord,
                    coord_vec, shapeFunctionType, centre, Dmat, B2);

                auto [e0, errorEstimates_e0, maxError_e0, minError_e0] =
                    matrix_integrator(-1.0, 1.0, f_bind_e0, rows, cols);

                auto [e1, errorEstimates_e1, maxError_e1, minError_e1] =
                    matrix_integrator(-1.0, 1.0, f_bind_e1, rows, cols);

                auto [e2, errorEstimates_e2, maxError_e2, minError_e2] =
                    matrix_integrator(-1.0, 1.0, f_bind_e2, rows, cols);

//#pragma omp critical
//                {
//                    std::cout << "e0: \n" << e0 << std::endl;
//                    std::cout << "e1: \n" << e1 << std::endl;
//                    std::cout << "e2: \n" << e2 << std::endl;
//                }
                e0_list_private.push_back(e0);
                e1_list_private.push_back(e1);
                e2_list_private.push_back(e2);
            }
            // Combine the private vectors into the shared vectors
#pragma omp critical
            e0_list.insert(e0_list.end(), e0_list_private.begin(),
                           e0_list_private.end());
#pragma omp critical
            e1_list.insert(e1_list.end(), e1_list_private.begin(),
                           e1_list_private.end());
#pragma omp critical
            e2_list.insert(e2_list.end(), e2_list_private.begin(),
                           e2_list_private.end());
        }


        auto stop = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time taken by calculation: " << duration.count()
                  << " microseconds" << std::endl;
        std::cout << "Time taken by calculation : "
                  << ((double)duration.count()) / 1000000 << " seconds"
                  << std::endl;

    }
}
