
#pragma once

#include "Filter.h"

#include <petsc.h>


class FilterWrapper
{
    public:
        // FilterWrapper(){}

        FilterWrapper(Filter *filter, PetscInt numConstraints)
        {
            this->filter         = filter;
            this->numConstraints = numConstraints;
        }

        virtual PetscErrorCode filterDesignVariable(Vec unfiltered, Vec filtered);
        virtual PetscErrorCode filterSensitivities(Vec position, Vec sensitivities, Vec* constraintSensitivities);

    private:
        Filter *filter;

        PetscInt numConstraints;

};
