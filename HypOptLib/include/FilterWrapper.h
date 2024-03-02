/***************************************************************************//**
 * @file FilterWrapper.h
 *
 * This file describes the layout of the FilterWrapper class.
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

#include <petsc.h>
#include "Filter.h"

/**
 * Class used to abstract filtering from the hyperoptimization class.
 *
 * Any inherited class can be used in lieu of this class when initializing hyperoptimization.
 * Instances are passed by reference, so any virtual functions in the derived class which
 * overwrite the ones provided here will be called by the design loop.
 */
class FilterWrapper
{
    public:
        /**
         * Empty constructor to allow passing by reference.
         *
         * @warning Do not instantiate using this constructor! This is only
         * provided as a means to pass by reference to the Hyperoptimization class.
         */
        FilterWrapper(){}

        /**
         * Correct constructor to use.
         *
         * @param filter Topopt filter object to be abstracted.
         * @param numConstraints Initialized by the TopOpt class: TopOpt::m.
         */
        FilterWrapper(Filter *filter, PetscInt numConstraints);

        /**
         * Applies the Topopt standard filter to the given design variables.
         *
         * @warning Unfiltered and filtered must be different objects (Petsc limitation). ie
         * cannot call FilterWrapper::filterDesignVariable(y, y).
         *
         * @param unfiltered the design variable to filter.
         * @param filtered [out] filtered positions. Must be identical settings to unfiltered Vec.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        virtual PetscErrorCode filterDesignVariable(Vec unfiltered, Vec filtered);

        /**
         * Applies the Topopt standard gradient filter to the provided sensitivities.
         *
         * @param position Position vector, should typically be filtered before using.
         * @param sensitivities [in/out] The sensitivities to filter. Function will overwrite this variable with the filtered result.
         * @param constraintSensitivities [in/out] The constraint sensitivities to filter. Function will overwrite this variable with
         * the filtered result. Should just be passed by reference, not an array.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        virtual PetscErrorCode filterSensitivities(Vec position, Vec sensitivities, Vec* constraintSensitivities);

    private:
        Filter *filter;
        PetscInt numConstraints;

};
