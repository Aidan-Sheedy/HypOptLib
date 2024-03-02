/***************************************************************************//**
 * @file SensitivitiesWrapper.h
 *
 * This file describes the layout of the SensitivitiesWrapper class.
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

#pragma once

#include "LinearElasticity.h"
#include "HypOptException.h"

#include <petsc.h>

/**
 * Class used to abstract sensitivities and objective functions from the hyperoptimization class.
 *
 * Any inherited class can be used in lieu of this class when initializing hyperoptimization.
 * Instances are passed by reference, so any virtual functions in the derived class which
 * overwrite the ones provided here will be called by the design loop.
 */
class SensitivitiesWrapper
{
    public:
        /**
         * Empty constructor to allow passing by reference.
         *
         * @warning Do not instantiate using this constructor! This is only
         * provided as a means to pass by reference to the Hyperoptimization class.
         */
        SensitivitiesWrapper(){}

        /**
         * Correct constructor to use.
         *
         * @param physics Topopt filter object to be abstracted.
         * @param Emin used in the LinearElasticity AssembleStiffnessMatrix function.
         * @param Emax used in the LinearElasticity AssembleStiffnessMatrix function.
         * @param penal optimization penalty power.
         * @param volfrac volume fraction.
         */
        SensitivitiesWrapper(LinearElasticity *physics,
                             PetscScalar Emin,
                             PetscScalar Emax,
                             PetscScalar penal,
                             PetscScalar volfrac)
        {
            if (nullptr == physics)
            {
                throw HypOptException("Null pointer passed to SensitivitiesWrapper.");
            }

            this->physics   = physics;
            this->Emin      = Emin;
            this->Emax      = Emax;
            this->penal     = penal;
            this->volfrac   = volfrac;
        }

        /**
         * Accessor for the Petsc solver iteration count.
         *
         * This is updated every time computeSensitivities or computeObjectiveFunction is called.
         *
         * @returns the number of iterations needed for the last finite element analysis calculation.
         */
        PetscInt getSolverIterationCount()
        {
            return petscSolverIterationCount;
        }

        /**
         * Wraps the LinearElasticity::ComputeSensitivities function.
         *
         * This solves the finite element analysis problem required to calcualte the sensitivities
         * and constraint sensitivities.
         *
         * @param filteredPositions filtered design particle position vector.
         * @param sensitivities [out] resulting sensitivities vector.
         * @param constraintSensitivities [out] resulting constraint sensitivities vector.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        virtual PetscErrorCode computeSensitivities(Vec filteredPositions,
                                                    Vec sensitivities,
                                                    Vec constraintSensitivities);

        /**
         * Wraps the LinearElasticity::ComputeObjectiveConstraints function.
         *
         * Note that the constraints are also calculated, but not used or returned. The finite element
         * analysis is calculated here too, which is the bottleneck in most hyperoptimization applications.
         *
         * @param filteredPositions filtered design particle position vector.
         * @param objectiveFunction [out] the resulting objective function.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        virtual PetscErrorCode computeObjectiveFunction(Vec filteredPositions, PetscScalar *objectiveFunction);
    private:
        LinearElasticity *physics;
        PetscScalar Emin;
        PetscScalar Emax;
        PetscScalar penal;
        PetscScalar volfrac;
        PetscInt    petscSolverIterationCount = 0;
};
