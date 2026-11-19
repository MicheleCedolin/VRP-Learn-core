cmake_minimum_required(VERSION 3.16)
project(vrplearn VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# -ffp-contract=off: no fused multiply-add, so every floating-point result is the plain IEEE-754
# double result.  This is what makes runs bit-for-bit reproducible across compilers and CPUs
# (and identical to the reference results the port was verified against).
# Do not add -ffast-math or -march=native without re-running the reference tests.
set(VRP_FLAGS -Wall -Wextra -ffp-contract=off)

add_library(vrplearn STATIC
  src/instance.cpp
  src/solution.cpp
  src/construct.cpp
  src/localsearch.cpp
  src/metaheuristics.cpp
  src/methods.cpp
  src/exact.cpp
  src/benchmarks.cpp)
target_include_directories(vrplearn PUBLIC include)
target_compile_options(vrplearn PRIVATE ${VRP_FLAGS})

add_executable(vrp apps/vrp.cpp)
target_link_libraries(vrp PRIVATE vrplearn)
target_compile_options(vrp PRIVATE ${VRP_FLAGS})

add_executable(vrp_experiments apps/experiments.cpp)
target_link_libraries(vrp_experiments PRIVATE vrplearn)
target_compile_options(vrp_experiments PRIVATE ${VRP_FLAGS})

enable_testing()
add_subdirectory(tests)
