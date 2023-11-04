
#pragma once

#include "LinearElasticity.h"
#include "TopOpt.h"

#include <petsc.h>

class SensitivitiesWrapper
{
    public:
        SensitivitiesWrapper(){}

        SensitivitiesWrapper(LinearElasticity *physics, TopOpt *opt)
        {
            this->physics = physics;
            this->opt     = opt;
        }

        virtual PetscErrorCode computeSensitivities(Vec filteredPositions,
                                                    Vec sensitivities,
                                                    Vec constraintSensitivities);
        virtual PetscErrorCode computeObjectiveFunction(Vec positions, PetscScalar *objectiveFunction);
    private:
        LinearElasticity *physics;
        TopOpt *opt;

};
