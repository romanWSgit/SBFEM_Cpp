# Add test executable
add_executable(SBFEM_Tests test_main.cpp test_math.cpp)


include_directories(${PROJECT_SOURCE_DIR}/include)


# Link the test executable against the Boost Test library and your source library
target_link_libraries(SBFEM_Tests
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        SBFEM_SRC
)
