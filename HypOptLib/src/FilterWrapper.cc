/***************************************************************************//**
 * @file FilterWrapper.cc
 *
 * This file implements the FilterWrapper class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#include "FilterWrapper.h"
#include "HypOptException.h"

FilterWrapper::FilterWrapper(Filter *filter, PetscInt numConstraints)
{
    if (nullptr == filter)
    {
        throw HypOptException("Invalid parameter: filter is null.");
    }

    this->filter         = filter;
    this->numConstraints = numConstraints;
}

PetscErrorCode FilterWrapper::filterDesignVariable(Vec unfiltered, Vec filtered)
{
    if (nullptr == filter)
    {
        throw HypOptException("Invalid Filter Wrapper! FilterWrapper::filter is null. You must initialize the FilterWrapper with the non-default constructor.");
    }

    Vec unused;
    PetscCall(VecDuplicate(unfiltered, &unused));
    filter->FilterProject(unfiltered, filtered, unused, PETSC_FALSE, 0, 0);
    PetscCall(VecDestroy(&unused));

    return 0;
}

PetscErrorCode FilterWrapper::filterSensitivities(Vec position, Vec sensitivities, Vec* constraintSensitivities)
{
    if (nullptr == filter)
    {
        throw HypOptException("Invalid Filter Wrapper! FilterWrapper::filter is null. You must initialize the FilterWrapper with the non-default constructor.");
    }

    Vec unused;
    PetscCall(VecDuplicate(position, &unused));
    PetscCall(filter->Gradients(position,
                                unused,
                                sensitivities,
                                numConstraints,
                                constraintSensitivities,
                                PETSC_FALSE,
                                0,
                                0));

    PetscCall(VecDestroy(&unused));

    return 0;
}
