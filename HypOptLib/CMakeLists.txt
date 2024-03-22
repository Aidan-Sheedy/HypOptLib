# Author: Thijs Smit, Dec 2020
# Updated: Aidan Sheedy, February 2024
# Copyright (C) 2020-2024 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

# check minimum required cmake version
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Command line options
option(InstallPythonSysPath "Will use the Python_SITELIB as the libray install directory." OFF)
option(CustomInstallPath "Sets a custom install path for the library." OFF)

# set the compiler
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# set the project name
project(HypOptLib VERSION 0.0.1
                  DESCRIPTION "Hyperoptimization extension of the functionality of the TopOpt_in_PETSc framework")

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)

# find_package(PythonInterp 3.7 REQUIRED)
# find_package(PythonLibs 3.7 REQUIRED)

## Set these values to point directly to your petsc location, if it is not found automatically.
# get_filename_component(PETSC_ABS_DIR "../../petsc" ABSOLUTE)
# set(PETSC_DIR /usr/local/lib/petsc)
# set(PETSC_ARCH "arch-linux-c-debug")

find_package(pybind11 REQUIRED)
find_package(PETSc REQUIRED)

message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

# Set install settings
set(CMAKE_INSTALL_RPATH "${PETSC_DIR}/${PETSC_ARCH}/lib")

add_definitions (-std=c++0x)

# include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_BINARY_DIR}/include" "${PETSC_INCLUDES}" "pybind11/include")

add_definitions (${PETSC_DEFINITIONS})
add_definitions (-g)

message (STATUS "PETSC_DEFINITIONS ${PETSC_DEFINITIONS}")

file(GLOB SOURCES "src/*.cc")

# Pybind11 library instructions
pybind11_add_module(HypOptLib SHARED ${SOURCES})
target_include_directories(HypOptLib PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include ${PETSC_INCLUDES} pybind11/include)
target_link_libraries(HypOptLib PUBLIC ${PETSC_LIBRARIES})

set_target_properties(
    HypOptLib
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "HypOptLib"
        LIBRARY_OUTPUT_DIRECTORY "../../run"
    )

# Setup install path based on command line options.
if(${InstallPythonSysPath})
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
    install(TARGETS HypOptLib DESTINATION ${Python_SITELIB})
    message (STATUS "InstallPythonSysPath is TRUE")
    message (STATUS "Installing to ${CMAKE_INSTALL_PREFIX}/${Python_SITELIB}")
elseif (NOT ("OFF" MATCHES "${CustomInstallPath}"))
    INSTALL(TARGETS HypOptLib DESTINATION ${CustomInstallPath})
    message (STATUS "Custom install path provided.")
    message (STATUS "Installing to ${CMAKE_INSTALL_PREFIX}/${CustomInstallPath}")
else()
    INSTALL(TARGETS HypOptLib DESTINATION lib)
endif()
