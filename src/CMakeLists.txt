cmake_minimum_required(VERSION 3.21)
project(sympleq LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CUDA_TOOLKIT_VERSION "12.4" CACHE STRING "CUDA Toolkit version")
option(CUDA_BUILD "Enable CUDA support" OFF)

message(STATUS "####### BUILD INFO #######")
message(STATUS "CMake preset: $ENV{CMAKE_GENERATOR_PRESET}")
message(STATUS "Adding: ${LIGHT_SRC_DIRECTORY}")
message(STATUS "Building for: ${CMAKE_BUILD_TYPE}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Cmake runtime output dir: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Enable vcpkg manifest mode
if(DEFINED ENV{VCPKG_ROOT})
    message(STATUS "Using vcpkg from $ENV{VCPKG_ROOT}")
endif()

set(PACKAGE_PY_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} CACHE STRING "Python package relative install path")
set(LIBS_PATH "${CMAKE_SOURCE_DIR}/sympleq.libs" CACHE STRING "Shared libs location")
set(PYTHON_PACKAGE_INSTALL_DIR "${CMAKE_BINARY_DIR}/install/${PACKAGE_PY_PATH}")

# Direct global output into the staged install tree for Release config only
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PYTHON_PACKAGE_INSTALL_DIR}")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PYTHON_PACKAGE_INSTALL_DIR}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${PYTHON_PACKAGE_INSTALL_DIR}")
endif()

# CUDA
if(CUDA_BUILD)
    message(STATUS "CUDA Toolkit version: ${CUDA_TOOLKIT_VERSION}")
    find_package(CUDAToolkit ${CUDA_TOOLKIT_VERSION} REQUIRED)
endif()

# Find GTest via vcpkg
find_package(GTest REQUIRED)

add_subdirectory("${LIGHT_SRC_DIRECTORY}")
