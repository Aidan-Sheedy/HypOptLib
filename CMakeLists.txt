# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

# check minimum required cmake version
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the compiler
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# set the project name
project(hypopt VERSION 0.0.1
                  DESCRIPTION "Hyperoptimization extension of the functionality of the TopOpt_in_PETSc framework")

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)

# find_package(PythonInterp 3.7 REQUIRED)
# find_package(PythonLibs 3.7 REQUIRED)

set(PETSC_DIR "/mnt/g/Projects/ENPH455/PETSc_paper_github/petsc")
set(PETSC_ARCH "arch-linux-c-opt")

find_package(PETSc REQUIRED)

message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

add_definitions (-std=c++0x)
include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_BINARY_DIR}/include" "${PETSC_INCLUDES}")

add_definitions (${PETSC_DEFINITIONS})
add_definitions (-g)

message (STATUS "PETSC_DEFINITIONS ${PETSC_DEFINITIONS}")

file(GLOB SOURCES "src/*.cc")
add_executable(hypopt ${SOURCES})

target_link_libraries(hypopt ${PETSC_LIBRARIES})
# target_link_libraries(topoptlib ${PYTHON_LIBRARIES})

set_target_properties(
    hypopt
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "hypopt"
        RUNTIME_OUTPUT_DIRECTORY "../bin"
    )
