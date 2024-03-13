/***************************************************************************//**
 * @file SensitivitiesWrapper.cc
 *
 * This file implements the default Sensitivities wrapper.
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

#include "SensitivitiesWrapper.h"

#include "PetscExtensions.h"

PetscErrorCode SensitivitiesWrapper::computeSensitivities(Vec filteredPositions,
                                                          Vec sensitivities,
                                                          Vec constraintSensitivities)
{
    PetscCall(
        physics->ComputeSensitivities(sensitivities,
                                      constraintSensitivities,
                                      filteredPositions,
                                      Emin,
                                      Emax,
                                      penal,
                                      volfrac,
                                      &petscSolverIterationCount);
    );

    return 0;
}

PetscErrorCode SensitivitiesWrapper::computeObjectiveFunction(Vec filteredPositions, PetscScalar *objectiveFunction)
{
    PetscScalar unused;
    PetscCall(
        physics->ComputeObjectiveConstraints(objectiveFunction,
                                            &unused,
                                            filteredPositions,
                                            Emin,
                                            Emax,
                                            penal,
                                            volfrac,
                                            &petscSolverIterationCount);
    );

    return 0;
}
