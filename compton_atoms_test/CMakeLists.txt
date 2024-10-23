cmake_minimum_required(VERSION 3.12)
project(full_test VERSION 1.0.0)

# Set some paths for the compiler
set(Boost_INCLUDE_DIR /Users/yashgursahani/Desktop/software/boost_1_82_0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin) # where to put the executable
set(OUTPUT_DIR /Users/yashgursahani/Desktop/UMD_Research/high_z_agn/mc_sim_output/logE_optimization_tests/) # where to put the output files

# Set the optimization level and C++ standard
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2") # set any flags, including optimization
set(CMAKE_BUILD_TYPE Release) # set the build type to release
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg") # set any flags for the linker, in this case a profiler

set(CMAKE_CXX_STANDARD 14) # want to use C++14
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include the Boost library
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Add your source files here
add_executable(full_test full_test/full_test/full_test.cpp full_test/full_test/full_test_functions.cpp) # note that we DO NOT need to add the header file here! Only the source files (.cpp)

# Link against the Boost libraries
target_link_libraries(full_test ${Boost_LIBRARIES})

# Run your program with an environment variable for the output directory
add_custom_target(run_my_program_with_env ALL
    COMMAND ${CMAKE_COMMAND} -E env OUTPUT_DIR=${OUTPUT_DIR} full_test
    DEPENDS full_test
    COMMENT "Running full_test with environment variable to generate output"
)