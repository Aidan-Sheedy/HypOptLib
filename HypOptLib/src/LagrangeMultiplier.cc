/***************************************************************************//**
 * @file LagrangeMultiplier.cc
 *
 * This file implements the default LagrangeMultiplier wrapper.
 *
 * @author Aidan Sheedy
 *
 * Copyright (C) 2024 Aidan Sheedy
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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