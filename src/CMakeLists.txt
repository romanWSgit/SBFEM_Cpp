enable_language(CXX)
enable_language(Fortran)

# Add your source files
add_library(SBFEM_SRC sbfem_math.cpp SuperElementJson.cpp
        SuperElement.cpp
        helper_functions.cpp helper_functions.cpp MaterialData.cpp
        sbfem_functions.cpp plot.cpp reorder_schur.cpp
        z_mat.cpp GraphicsController.cpp NURBSVizalisation.cpp
        ShapeFunctionType.cpp
        sbfem_driver.cpp
        FortranInteroperability.cpp
)

# Add the Fortran source files
add_library(fortran_lib
        STATIC
        fortran/sbfem_functions.f90
        fortran/helper_functions.f90
        fortran/math.f90
        fortran/functional.f90)

# Include the header files
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external/tinynurbs/include)

target_link_libraries(SBFEM_SRC
        PRIVATE
        fast_matrix_market::fast_matrix_market
        nlohmann_json::nlohmann_json
        ${Boost_LIBRARIES}
        fortran_lib
        gfortran
        )

