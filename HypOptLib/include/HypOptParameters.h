/***************************************************************************//**
 * @file HypOptParameters.h
 *
 * This file describes structures useful for hyperoptimization.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#pragma once

#include <petsc.h>
#include <vector>
#include <set>

/**
 * Structure defining all aspects of the design space.
 *
 * Mostly a convenience structure, as these parameters are often passed together.
 */
typedef struct
{
    Vec position;
    Vec velocity;
    Vec evenNoseHooverPosition;
    Vec evenNoseHooverVelocity;
    Vec oddNoseHooverPosition;
    Vec oddNoseHooverVelocity;
} HypOptParameters;

/** 
 * Types of boundary conditions currently supported.
 */
typedef enum
{
    FIXED_POINT,
    LOAD,
} BoundaryConditionType;

/** 
 * Domain definition. These coordinates define the size of the problem.
 */
typedef struct
{
    PetscScalar xMinimum;
    PetscScalar xMaximum;
    PetscScalar yMinimum;
    PetscScalar yMaximum;
    PetscScalar zMinimum;
    PetscScalar zMaximum;
} DomainCoordinates;

/**
 * Boundary condition set by user. The types supported are fixed points and loads.
 */
typedef struct
{
    BoundaryConditionType    type;
    std::vector<PetscScalar> xRange;
    std::vector<PetscScalar> yRange;
    std::vector<PetscScalar> zRange;
    std::set<PetscInt>       degreesOfFreedom;
    PetscScalar              value; /** Only used for load types currently, typically around -0.001. */
} BoundaryCondition;

/**
 * Debug log verbosity levels.
 */
enum verbosity
{
    QUIET = 0,
    INFO,
    DEBUG,
};