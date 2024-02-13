//
// Created by Roman Wallner- Silberhuber on 02.06.23.
//

#ifndef SBFEM_MATERIALDATA_H
#define SBFEM_MATERIALDATA_H

#include <iostream>
#include <fstream>
#include <string>
#include <optional>
#include <nlohmann/json.hpp>

struct Params {
    double Em;
    double nu;
    double rho;
    double Ex;
    double Ey;
    double G;
    double omega;
    double c11;
    double c13;
    double c33;
    double c44;
    double e31;
    double e33;
    double e15;
    double epsilon11;
    double epsilon33;
};

struct MaterialData {
    std::string Material;
    std::string MaterialProp;
    std::string D;
    Params params;
};

struct StructData {
    std::string structType;
    std::vector<MaterialData> MaterialData;
};



Params to_params(const nlohmann::json& j);

MaterialData to_material_data(const nlohmann::json& j);

StructData to_struct_data(const nlohmann::json& j);

/**
 * \relates MaterialData
 * \fn     std::optional<MaterialData> findStructDataAndMaterial(const std::vector<StructData>& structs, \
                                                  const std::string& structType, const std::string& material);
 * \brief  useful to to access a specific material
 * \param  const std::vector<StructData>& structs Description of the first parameter
 * \param  const std::string& structType
 * \param  const std::string& material
 * \return returns an std::optional<MaterialData> object iff the field name is found in the json handle
 *
 * one can access to object by the "->" operator: the structure is: material_data: -> |
 *                                                                                    | - params
 *                                                                                    | - D
 *                                                                                    | - Material
 *                                                                                    | - MaterialProp
 */
std::optional<MaterialData> findStructDataAndMaterial(const std::vector<StructData>& structs, \
                                                         const std::string& structType, const std::string& material);

#endif //SBFEM_MATERIALDATA_H
