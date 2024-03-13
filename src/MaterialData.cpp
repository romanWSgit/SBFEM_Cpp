//
// Created by Roman Wallner- Silberhuber on 02.06.23.
//

#include "MaterialData.h"

// Function to convert a json object to a Params object
Params to_params(const nlohmann::json &j)
{
    Params params;
    params.Em = j.value("Em", 0.0);
    params.nu = j.value("ν", 0.0);
    params.rho = j.value("ρ", 0.0);
    params.Ex = j.value("Ex", 0.0);
    params.Ey = j.value("Ey", 0.0);
    params.G = j.value("G", 0.0);
    params.omega = j.value("ω", 0.0);
    params.c11 = j.value("c11", 0.0);
    params.c33 = j.value("c33", 0.0);
    params.c44 = j.value("c44", 0.0);
    params.e31 = j.value("e31", 0.0);
    params.e33 = j.value("e33", 0.0);
    params.e15 = j.value("e15", 0.0);
    params.epsilon11 = j.value("ϵ11", 0.0);
    params.epsilon33 = j.value("ϵ33", 0.0);

    return params;
}

// Function to convert a json object to a MaterialData object
MaterialData to_material_data(const nlohmann::json &j)
{
    MaterialData material_data;
    material_data.Material = j.value("Material", "");
    material_data.MaterialProp = j.value("MaterialProp", "");
    material_data.D = j.value("D", "");
    material_data.params = to_params(j["params"]);
    return material_data;
}

// Function to convert a json object to a StructData object
StructData to_struct_data(const nlohmann::json &j)
{
    StructData struct_data;
    struct_data.structType = j.value("struct", "");
    for (const auto &material_data_json : j["MaterialData"])
    {
        struct_data.MaterialData.push_back(
            to_material_data(material_data_json));
    }
    return struct_data;
}

std::optional<MaterialData> findStructDataAndMaterial(
    const std::vector<StructData> &structs, const std::string &structType,
    const std::string &material)
{
    for (const auto &s : structs)
    {
        if (s.structType == structType)
        {
            for (const auto &m : s.MaterialData)
            {
                if (m.Material == material)
                {
                    return m;
                }
            }
        }
    }
    return std::nullopt;
}


double replaceExpressions(const std::string &exprStr, const Params &params)
{
    double Em = params.Em; // Temporary variables
    double nu = params.nu;

    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_variable("Em", Em);
    symbol_table.add_variable("nu", nu);

    exprtk::expression<double> expression;
    expression.register_symbol_table(symbol_table);

    exprtk::parser<double> parser;
    if (!parser.compile(exprStr, expression))
    {
        throw std::runtime_error("Failed to parse expression!");
    }

    return expression.value();
}


Eigen::MatrixXd MaterialData::getDMatrix() const {
    std::string processedD = std::regex_replace(D, std::regex("ν"), "nu"); // Assuming regex issue is resolved
    std::istringstream iss(processedD);
    std::string token;
    std::vector<std::vector<double>> rows;

    // Skip initial {{
    std::getline(iss, token, '{');

    while (std::getline(iss, token, '{')) {
        std::vector<double> row;
        std::string cell;

        // Get the content until the closing '}' of the current row
        std::getline(iss, cell, '}');

        // Now, split the cell content by commas
        std::istringstream cellStream(cell);
        std::string value;

        while (std::getline(cellStream, value, ',')) {
            // Check if the value is empty or just a closing brace which should not happen now
            if (!value.empty() && value.find('}') == std::string::npos) {
                try {
                    double cellValue = replaceExpressions(value, params);
                    row.push_back(cellValue);
                } catch (const std::runtime_error& e) {
                    std::cerr << "Error parsing expression: " << value << ". Error: " << e.what() << std::endl;
                    throw; // Rethrow or handle as needed
                }
            }
        }

        if (!row.empty()) {
            rows.push_back(row);
        }
    }

    Eigen::MatrixXd mat(rows.size(), rows.empty() ? 0 : rows.front().size());
    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < rows[i].size(); ++j) {
            mat(i, j) = rows[i][j];
        }
    }

    return mat;
}