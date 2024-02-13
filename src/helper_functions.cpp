//
// Created by Roman Wallner- Silberhuber on 25.05.23.
//

#include "helper_functions.h"





Eigen::MatrixXd readMatrixFromFile (const std::string &filePath)
{
    Eigen::MatrixXd matrix;
    try
        {
            std::ifstream file (filePath);
            fast_matrix_market::read_matrix_market_eigen_dense (file, matrix);
            file.close ();
        }
    catch (const fast_matrix_market::invalid_mm &e)
        {
            std::cerr << "Error: " << e.what () << std::endl;
            // Handle the error appropriately
            // For example, you might want to exit the function or return an empty matrix
        }
    return matrix;
}

StructData getMaterialDataForType (std::string &search_materialType, \
                                      std::vector<StructData> &material_data_list)
{

    auto it = std::find_if (material_data_list.begin (), material_data_list.end (),
                            [&search_materialType] (const StructData &s)
                            {
                                return s.structType == search_materialType;
                            });

    if (it != material_data_list.end ())
        {
            // Found it
            StructData &found_material_data = *it;
            // Now you can use found_struct...
            std::cout << "Found: " << found_material_data.structType
                      << std::endl;
            return found_material_data;
        }
    else
        {
            // Didn't find it
            std::cout << "Did not find struct with type: "
                      << search_materialType << std::endl;
            throw std::runtime_error ("Material not found in dataset");
        }
}



