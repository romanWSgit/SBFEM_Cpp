
# Add your source files
add_library(SBFEM_SRC sbfem_math.cpp SuperElement.cpp
        helper_functions.cpp helper_functions.cpp MaterialData.cpp
        sbfem_functions.cpp plot.cpp reorder_schur.cpp
        z_mat.cpp)
# Include the header files
include_directories(${PROJECT_SOURCE_DIR}/include)

target_link_libraries(SBFEM_SRC
        PRIVATE
        fast_matrix_market::fast_matrix_market
        nlohmann_json::nlohmann_json
        ${Boost_LIBRARIES}
        )

