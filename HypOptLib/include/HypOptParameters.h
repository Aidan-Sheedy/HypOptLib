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
    FIXED_POINT,    /**< Cells are fixed at one. Future support may allow for arbitrary fixed point values and voids. */
    LOAD,           /**< Cells have a load applied at this ponit. */
} BoundaryConditionType;

/** 
 * Domain definition. These coordinates define the size of the problem.
 *
 * This allows for varying grid sizes without modifying the boundary conditions.
 * The domain is used in all boundary condition calculations.
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
 *
 * The range in each dimension are defined as (min, max), where max >= min. If max = min, then
 * the range is singular in that axis. If all ranges are singular, the boundary condition describes
 * a single piont.
 */
typedef struct
{
    BoundaryConditionType    type;              /**< fixed point or laod. */
    std::vector<PetscScalar> xRange;            /**< Range at which the boundary condition applies in the x direction. */
    std::vector<PetscScalar> yRange;            /**< Range at which the boundary condition applies in the y direction. */
    std::vector<PetscScalar> zRange;            /**< Range at which the boundary condition applies in the z direction. */
    std::set<PetscInt>       degreesOfFreedom;  /**< set of desired degrees of freedom. Must one or multiple of 0, 1, or 2,
                                                     corresponding to x, y, and z respectively. */
    PetscScalar              value;             /**< Only used for load types currently, typically around -0.001. If there is more
                                                     than one degree of freedom, the value is applies to each one respectively. */
} BoundaryCondition;

/**
 * Debug log verbosity levels.
 */
enum verbosity
{
    QUIET = 0,  /**< Most logs are muted. */
    INFO,       /**< Default log level. */
    DEBUG,      /**< Higher level logs are printed. */
};