/***************************************************************************//**
 * @file FilterWrapper.h
 *
 * This file describes the layout of the FilterWrapper class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#pragma once

#include <petsc.h>
#include "Filter.h"

/**
 * Class used to abstract filtering from the hyperoptimization class.
 *
 * Any inherited class can be used in lue of this class when initializing hyperoptimization.
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
         * provided as a means to pass by reference to the hyperoptimization class.
         *
         * @note Inherited classes are free to use this constructor as they may not need
         * to use the default parameters used here.
         */
        FilterWrapper(){}

        /**
         * Correct constructor to use.
         *
         * @param filter Topopt filter opbject to be abstracted.
         * @param numConstraints Initialized by the TopOpt class: TopOpt::m.
         */
        FilterWrapper(Filter *filter, PetscInt numConstraints);

        /**
         * Applies the Topopt standard filter to the given design variables.
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
