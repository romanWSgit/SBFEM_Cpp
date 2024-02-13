//
// Created by Roman Wallner- Silberhuber on 02.06.23.
//

#include "MaterialData.h"




// Function to convert a json object to a Params object
Params to_params(const nlohmann::json& j) {
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
MaterialData to_material_data(const nlohmann::json& j)
{
    MaterialData material_data;
    material_data.Material = j.value("Material", "");
    material_data.MaterialProp = j.value("MaterialProp", "");
    material_data.D = j.value("D", "");
    material_data.params = to_params(j["params"]);
    return material_data;
}

// Function to convert a json object to a StructData object
StructData to_struct_data(const nlohmann::json& j)
{
    StructData struct_data;
    struct_data.structType = j.value("struct", "");
    for (const auto& material_data_json : j["MaterialData"]) {
        struct_data.MaterialData.push_back(to_material_data(material_data_json));
    }
    return struct_data;
}


std::optional<MaterialData> findStructDataAndMaterial(const std::vector<StructData>& structs, \
                                                         const std::string& structType, const std::string& material)
{
    for (const auto& s : structs) {
        if (s.structType == structType) {
            for (const auto& m : s.MaterialData) {
                if (m.Material == material) {
                    return m;
                }
            }
        }
    }
    return std::nullopt;
}