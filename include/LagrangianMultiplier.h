
#pragma once

#include "TopOpt.h"
#include "Filter.h"

#include <petsc.h>

class LagrangianMultiplier
{
    public:
        LagrangianMultiplier(){}

        LagrangianMultiplier(Filter *filter, TopOpt *opt)
        {
            this->filter    = filter;
            this->opt       = opt;
        }

        PetscScalar computeLagrangianMultiplier(Vec positions, Vec C, PetscInt numParticles, PetscScalar *returnValue);

    private:
        Filter *filter;

        TopOpt *opt;

};