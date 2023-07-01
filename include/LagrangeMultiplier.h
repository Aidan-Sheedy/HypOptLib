
#pragma once

#include "TopOpt.h"
#include "Filter.h"

#include <petsc.h>

class LagrangeMultiplier
{
    public:
        LagrangeMultiplier(){}

        LagrangeMultiplier(Filter *filter, TopOpt *opt)
        {
            this->filter    = filter;
            this->opt       = opt;
        }

        virtual PetscScalar computeLagrangeMultiplier(Vec positions, Vec C, PetscInt numParticles, PetscScalar *returnValue);

    private:
        Filter *filter;

        TopOpt *opt;

};
