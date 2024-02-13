cmake_minimum_required(VERSION 3.24)
project(SBFEM)

set(CMAKE_CXX_STANDARD 17)

### set python Version:
#set(Python_INCLUDE_DIRS "/Users/roman_w_s/mambaforge/include")
#set(Python_LIBRARIES "/Users/roman_w_s/mambaforge/lib")
#set(Python_EXECUTABLE "/Users/roman_w_s/mambaforge/bin/python3")
#find_package(Python COMPONENTS Interpreter Development)



# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

#include(FetchContent)
#include(FetchContent)
#FetchContent_Declare(
#  googletest
#  # Specify the commit you depend on and update it regularly.
#  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
#)


find_package(Eigen3 REQUIRED)
if(NOT EXISTS ${EIGEN3_INCLUDE_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})
message("Eigen include dir: " ${EIGEN3_INCLUDE_DIR})

set(ENV{PATH} "/opt/local/bin:$ENV{PATH}")
find_program(GNUPLOT_EXECUTABLE NAMES gnuplot)
message(STATUS "Gnuplot executable: ${GNUPLOT_EXECUTABLE}")




find_package(nlohmann_json REQUIRED)

find_package(Matplot++ REQUIRED)

# Find the Boost library
find_package(Boost REQUIRED COMPONENTS
        unit_test_framework
        filesystem
        system
        iostreams )


include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/gnuplot-iostream")

# Include Boost in the project
include_directories(${Boost_INCLUDE_DIRS})

#find_package(pybind11 CONFIG REQUIRED)

## Find PLplot package
#find_package(PkgConfig REQUIRED)
#pkg_search_module(PLPLOT REQUIRED plplot-c++)
#include_directories(${PLPLOT_INCLUDE_DIRS})


include(FetchContent)
FetchContent_Declare(
        fast_matrix_market
        GIT_REPOSITORY https://github.com/alugowski/fast_matrix_market
        GIT_TAG main
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(fast_matrix_market)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    message("-- Found LAPACK library")
endif()

find_package(BLAS)
if(LAPACK_FOUND)
    message("-- Found BLAS library")
endif()

set(SLICOT_LIBRARY "/usr/local/lib/libslicot.a")


add_executable(SBFEM
        main.cpp
        include/reorder_schur.h
        include/exceptions.h
        include/z_mat.h)

include_directories(${PROJECT_SOURCE_DIR}/include)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

# Include the header files


target_link_libraries(SBFEM
        PRIVATE
#        ${VTK_LIBRARIES}
#        ${PLPLOT_LIBRARIES}
#        ${PLPLOT_LIBRARIES}
        fast_matrix_market::fast_matrix_market
        nlohmann_json::nlohmann_json
        ${Boost_LIBRARIES}
        SBFEM_SRC
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${SLICOT_LIBRARY}
        Matplot++::matplot
        )

# Define BOOST_TEST_DYN_LINK to enable dynamic linking
#target_compile_definitions(SBFEM PRIVATE "BOOST_TEST_DYN_LINK")

#.-------------------------------------------------------

# DOXYGEN
message("cmake dir :${CMAKE_CURRENT_SOURCE_DIR}")
# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
    message("Doxygen executable = ${DOXYGEN_EXECUTABLE}")
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

# Add tests
