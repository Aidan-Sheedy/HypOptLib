/***************************************************************************//**
 * @file SensitivitiesWrapper.h
 *
 * This file describes the layout of the SensitivitiesWrapper class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
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
};
