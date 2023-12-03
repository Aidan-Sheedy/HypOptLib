/***************************************************************************//**
 * @file LagrangeMultiplier.cc
 *
 * This file implements the default LagrangeMultiplier wrapper.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#include "LagrangeMultiplier.h"

PetscScalar LagrangeMultiplier::computeLagrangeMultiplier(Vec positions, Vec C, PetscInt numParticles, PetscScalar *LagrangeMultiplier)
{
    PetscErrorCode errorStatus = 0;
    PetscScalar positionSum;
    PetscScalar CSum;
    Vec filteredPosition;
    Vec filteredC;
    
    PetscCall(VecDuplicate(positions, &filteredPosition));
    PetscCall(VecDuplicate(C, &filteredC));

    /* Filter positions and C */
    PetscCall(this->filter.filterDesignVariable(positions, filteredPosition));
    PetscCall(this->filter.filterDesignVariable(C, filteredC));

    /* Calculate sums */
    PetscCall(VecSum(filteredPosition, &positionSum));
    PetscCall(VecSum(filteredC, &CSum));

    /* Calculate Lagrangee multiplier */
    *LagrangeMultiplier = (numParticles * volfrac - positionSum )/ CSum;

    PetscCall(VecDestroy(&filteredPosition));
    PetscCall(VecDestroy(&filteredC));

    return errorStatus;
}