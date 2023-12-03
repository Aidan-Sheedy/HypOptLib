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
