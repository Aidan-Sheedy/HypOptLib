/***************************************************************************//**
 * @file LagrangeMultiplier.h
 *
 * This file describes the layout of the LagrangeMultiplier wrapper.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#pragma once

#include "FilterWrapper.h"

#include <petsc.h>

class LagrangeMultiplier
{
    public:
        /**
         * Empty constructor to allow passing by reference.
         *
         * @warning Do not instantiate using this constructor! This is only
         * provided as a means to pass by reference to the Hyperoptimization class.
         */
        LagrangeMultiplier(){}

        /**
         * Correct constructor to use.
         *
         * @param filter abstracted filter wrapper. Note it is passed by reference to allow for any arbitrary filter implementation.
         * @param volfrac target volume fraction of the system.
         */
        LagrangeMultiplier(FilterWrapper& filter, PetscScalar volfrac)
        {
            this->filter    = filter;
            this->volfrac   = volfrac;
        }

        /**
         * Default Lagrange multiplier calculation.
         *
         * For any given problem, this function can be overwritten by derived classes to implement the Lagrange Multiplier calculation
         * specific to said problem. This particular function implements the following equation:
         *
         * @f[
         * \lambda = \frac{N c_v - \sum_i^N{\tilde{x_i}}}{\sum_i^N{\tilde{C_i}}},
         * @f]
         *
         * where:
         *  - \f$ \lambda \f$ is the Lagrange multiplier,
         *  - \f$ N \f$ is the number of design particles in the system,
         *  - \f$ c_v \f$ is volume fraction,
         *  - \f$ \tilde{x_i} \f$ is the i'th filtered design variable position,
         *  - \f$ \tilde{C_i} \f$ is the i'th filtered C vector element.
         *
         * @todo See the hyperoptimization summary for the C vector definition.
         *
         * @param positions vector of design variable positions.
         * @param C the C vector.
         * @param numParticles the number of design variable particles.
         * @param LagrangeMultiplier [out] the resulting Lagrange Multiplier.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        virtual PetscScalar computeLagrangeMultiplier(Vec positions, Vec C, PetscInt numParticles, PetscScalar *LagrangeMultiplier);

    private:
        FilterWrapper filter;
        PetscScalar volfrac;
};
