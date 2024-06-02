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
option(CUSTOM_INSTALL_PREFIX "Sets a custom prefix path for the library installation." OFF)
option(FORCE_INSTALL_PATH "Sets a custom install path for the library." OFF)
option(USE_PETSC_RPATH "Sets the RPATH to use the Petsc library directory." OFF)
option(USE_PYTHON_VERSION "The install location will try to use this python version if it is found." OFF)


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

MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children LIST_DIRECTORIES true RELATIVE ${curdir} ${curdir}/*)
    SET(dirlist "")
    FOREACH(child ${children})
        # message(STATUS "    ${curdir} child: ${child}")
        if (IS_DIRECTORY ${curdir}/${child})
            # message(STATUS "       ${curdir}/${child} IS DIRECTORY!")
            LIST(APPEND dirlist ${child})
        else()
            # message(STATUS "       ${curdir}/${child} isn't :(")
        endif()
    ENDFOREACH()
    set(${result} ${dirlist})
    # message(STATUS "dirlist: ${dirlist}")
    # message(STATUS "result: ${result}")
ENDMACRO()

MACRO(FIND_INSTALLPATH return)
    set(INSTALLPATH /lib)
    if(EXISTS ${CMAKE_INSTALL_PREFIX}${INSTALLPATH})
        SUBDIRLIST(PREFIX_SUBDIRS ${CMAKE_INSTALL_PREFIX}${INSTALLPATH})
        foreach(subdir ${PREFIX_SUBDIRS})
            if (${subdir} MATCHES python.*)
                if (${subdir} MATCHES /^python$/)
                    set(INSTALLPATH /lib/${subdir})
                    break()
                elseif (EXISTS ${CMAKE_INSTALL_PREFIX}lib/${subdir}/site-packages)
                    set(INSTALLPATH /lib/${subdir}/site-packages)
                    break()
                endif()
            endif()
        endforeach()
    endif()
    set(${return} ${INSTALLPATH})
ENDMACRO()

# Setup install path based on command line options.
if (NOT ("OFF" MATCHES "${CUSTOM_INSTALL_PREFIX}"))
    set(CMAKE_INSTALL_PREFIX ${CUSTOM_INSTALL_PREFIX})
endif()


if (NOT ("OFF" MATCHES "${FORCE_INSTALL_PATH}"))
    set(INSTALLPATH ${FORCE_INSTALL_PATH})
else()
    if (NOT ("OFF" MATCHES "${USE_PYTHON_VERSION}"))
        if(EXISTS ${CMAKE_INSTALL_PREFIX}lib/${USE_PYTHON_VERSION}/site-packages)
            set(INSTALLPATH lib/${USE_PYTHON_VERSION}/site-packages)
        else()
            message(WARNING "Could not find python version: ${USE_PYTHON_VERSION}")
            FIND_INSTALLPATH(INSTALLPATH)
        endif()
    else()
        FIND_INSTALLPATH(INSTALLPATH)
    endif()

    if (NOT EXISTS ${CMAKE_INSTALL_PREFIX}/${INSTALLPATH})
        message(ERROR "Could not find appropriate install path in ${CMAKE_INSTALL_PREFIX}")
        return()
    endif()
endif()

message(STATUS "Install path: ${INSTALLPATH}")

# Set install settings
if ("ON" MATCHES "${USE_PETSC_RPATH}")
    set(CMAKE_INSTALL_RPATH "${PETSC_DIR}/${PETSC_ARCH}/lib")
else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}${INSTALLPATH}")
endif()

## This line may be needed if petsc is not being found or if the library is not being imported
## to python after building with make install.
# list( APPEND CMAKE_INSTALL_RPATH ${PETSC_DIR}/lib )

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

if(${InstallPythonSysPath})
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
    install(TARGETS HypOptLib DESTINATION ${Python_SITELIB})
    message (STATUS "InstallPythonSysPath is TRUE")
    message (STATUS "Installing to ${CMAKE_INSTALL_PREFIX}${Python_SITELIB}")
else()
    message (STATUS "installing to ${CMAKE_INSTALL_PREFIX}${INSTALLPATH}")
    message(STATUS "    Prefix path: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "    Install path: ${INSTALLPATH}")
    INSTALL(TARGETS HypOptLib DESTINATION ${CMAKE_INSTALL_PREFIX}${INSTALLPATH})
endif()
