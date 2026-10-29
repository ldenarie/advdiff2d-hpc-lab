cmake_minimum_required(VERSION 3.22)
project(advdiff2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ADVDIFF_BUILD_TESTS "Build the test suite" ON)
option(ADVDIFF_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

# Core library: data structures and kernels.
add_library(advdiff src/kernel_seq.cpp)
target_include_directories(advdiff PUBLIC src)
target_compile_options(advdiff PUBLIC
  -Wall -Wextra -Wpedantic
  $<$<BOOL:${ADVDIFF_WARNINGS_AS_ERRORS}>:-Werror>)

# Driver.
add_executable(advdiff2d src/main.cpp)
target_link_libraries(advdiff2d PRIVATE advdiff)

# Tests (Catch2 v3 from the system package catch2).
if(ADVDIFF_BUILD_TESTS)
  enable_testing()
  find_package(Catch2 3 REQUIRED)
  include(Catch)
  add_executable(test_kernel tests/test_kernel.cpp)
  target_link_libraries(test_kernel PRIVATE advdiff Catch2::Catch2WithMain)
  # Catch2 exits with code 4 when a test is skipped; report it as skipped, not failed.
  catch_discover_tests(test_kernel PROPERTIES SKIP_RETURN_CODE 4)
endif()
