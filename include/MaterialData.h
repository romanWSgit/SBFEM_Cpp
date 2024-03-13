//
// Created by Roman Wallner- Silberhuber on 02.06.23.
//

#ifndef SBFEM_MATERIALDATA_H
#define SBFEM_MATERIALDATA_H

#include "exprtk.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <optional>
#include <regex>
#include <string>

struct Params
{
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

struct MaterialData
{
    std::string Material;
    std::string MaterialProp;
    std::string D;
    Params params;

    /**
     * @brief Get the D matrix.
     *
     * This function parses the string D and converts it into a matrix using the
     * given parameters.
     *
     * @return The D matrix.
     * @throws std::runtime_error if there is an error parsing the expression.
     */
    Eigen::MatrixXd getDMatrix() const;
};

struct StructData
{
    std::string structType;
    std::vector<MaterialData> MaterialData;
};

Params to_params(const nlohmann::json &j);

MaterialData to_material_data(const nlohmann::json &j);

StructData to_struct_data(const nlohmann::json &j);

/**
 * \relates MaterialData
 * \fn     std::optional<MaterialData> findStructDataAndMaterial(const
 std::vector<StructData>& structs, \ const std::string& structType, const
 std::string& material);
 * \brief  useful to to access a specific material
 * \param  const std::vector<StructData>& structs Description of the first
 parameter
 * \param  const std::string& structType
 * \param  const std::string& material
 * \return returns an std::optional<MaterialData> object iff the field name is
 found in the json handle
 *
 * one can access to object by the "->" operator: the structure is:
 material_data: -> |
 * | - params
 * | - D
 * | - Material
 * | - MaterialProp
 */
std::optional<MaterialData> findStructDataAndMaterial(
    const std::vector<StructData> &structs, const std::string &structType,
    const std::string &material);

/**
 * @brief Replaces the expressions in the given string using the provided
 * parameters.
 *
 * This function replaces expressions enclosed within curly braces in the given
 * string with their corresponding values from the provided parameters. The
 * expressions must follow the format "{key}" where 'key' represents a valid
 * parameter key.
 *
 * @param exprStr The string containing expressions to be replaced.
 * @param params A key-value map of parameters.
 * @return The resulting string after replacing the expressions.
 *
 * @note This function does not modify the original string.
 */
double replaceExpressions(const std::string &exprStr, const Params &params);

#endif // SBFEM_MATERIALDATA_H
