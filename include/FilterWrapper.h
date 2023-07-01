
#pragma once

#include "Filter.h"

#include <petsc.h>


class FilterWrapper
{
    public:
        FilterWrapper(){}

        FilterWrapper(Filter *filter)
        {
            this->filter    = filter;
        }

        virtual PetscErrorCode filterDesignVariable(Vec unfiltered, Vec filtered);
        virtual PetscErrorCode filterSensitivities(Vec position, Vec sensitivities, Vec* constraintSensitivities);

    private:
        Filter *filter;

};
