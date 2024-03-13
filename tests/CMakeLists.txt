# Add test executable


add_executable(SBFEM_Tests
        test_main.cpp
        test_math.cpp
        test_sbfem_fct.cpp)


include_directories(${PROJECT_SOURCE_DIR}/include)


## Link the test executable against the Boost Test library and your source library
#target_link_libraries(SBFEM_Tests
#
#)

target_link_libraries(SBFEM_Tests PUBLIC
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        SBFEM_SRC
        fast_matrix_market::fast_matrix_market)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message("OpenMp found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    target_link_libraries(SBFEM_Tests PUBLIC OpenMP::OpenMP_CXX)
endif()