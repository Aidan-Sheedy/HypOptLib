/***************************************************************************//**
 * @file FilterWrapper.cc
 *
 * This file implements the FilterWrapper class.
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
