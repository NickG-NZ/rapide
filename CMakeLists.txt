cmake_minimum_required(VERSION 3.23)
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

option(RAPIDE_WITH_COAL "Build the Coal collision adapter and contact tests" ON)
if(RAPIDE_WITH_COAL)
    find_package(coal 3.0.4 REQUIRED CONFIG)
    add_library(rapide_collision INTERFACE)
    target_link_libraries(rapide_collision INTERFACE rapide coal::coal)
endif()

# Examples
add_executable(motorcycle_rollout examples/motorcycle_rollout.cpp)
target_link_libraries(motorcycle_rollout PRIVATE rapide)
target_compile_options(motorcycle_rollout PRIVATE -Wall -Wextra -O2)

add_executable(four_bar_rollout examples/four_bar_rollout.cpp)
target_link_libraries(four_bar_rollout PRIVATE rapide)
target_compile_options(four_bar_rollout PRIVATE -Wall -Wextra)

add_executable(wheel_rollout examples/wheel_rollout.cpp)
target_link_libraries(wheel_rollout PRIVATE rapide)
target_compile_options(wheel_rollout PRIVATE -Wall -Wextra)

add_executable(bike_rollout examples/bike_rollout.cpp)
target_link_libraries(bike_rollout PRIVATE rapide)
target_compile_options(bike_rollout PRIVATE -Wall -Wextra)

# Tests
include(CTest)
if(BUILD_TESTING)
    find_package(Catch2 3 REQUIRED CONFIG)
    foreach(test se3 spatial aba tire maximal planar_tire bike integration backends)
        add_executable(test_${test} test/test_${test}.cpp)
        target_link_libraries(test_${test} PRIVATE rapide Catch2::Catch2WithMain)
        target_compile_options(test_${test} PRIVATE -Wall -Wextra -O2)
        add_test(NAME ${test} COMMAND test_${test})
    endforeach()
    if(RAPIDE_WITH_COAL)
        add_executable(test_collision test/test_collision.cpp)
        target_link_libraries(test_collision PRIVATE rapide_collision Catch2::Catch2WithMain)
        add_test(NAME collision COMMAND test_collision)
    endif()
endif()
