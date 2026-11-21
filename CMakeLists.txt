cmake_minimum_required(VERSION 3.20)
project(rapide VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# rapide is effectively header-only (algorithms are templated).
# The library target exists for linking convenience and compile-flag propagation.
add_library(rapide INTERFACE)
target_include_directories(rapide INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(rapide INTERFACE Eigen3::Eigen)

# Examples
add_executable(motorcycle_rollout examples/motorcycle_rollout.cpp)
target_link_libraries(motorcycle_rollout PRIVATE rapide)
target_compile_options(motorcycle_rollout PRIVATE -Wall -Wextra -O2)

# Tests
find_package(Catch2 3 QUIET)
if(Catch2_FOUND)
    enable_testing()
    foreach(test se3 spatial aba tire)
        add_executable(test_${test} test/test_${test}.cpp)
        target_link_libraries(test_${test} PRIVATE rapide Catch2::Catch2WithMain)
        target_compile_options(test_${test} PRIVATE -Wall -Wextra -O2)
        add_test(NAME ${test} COMMAND test_${test})
    endforeach()
endif()
