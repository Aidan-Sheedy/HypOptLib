#include "SensitivitiesWrapper.h"

#include "PetscExtensions.h"

PetscErrorCode SensitivitiesWrapper::computeSensitivities(Vec filteredPositions,
                                                          Vec sensitivities,
                                                          Vec constraintSensitivities)
{
    PetscCall(
        physics->ComputeSensitivities(sensitivities,
                                      constraintSensitivities,
                                      filteredPositions,
                                      opt->Emin,
                                      opt->Emax,
                                      opt->penal,
                                      opt->volfrac);
    );

    return 0;
}

PetscErrorCode SensitivitiesWrapper::computeObjectiveFunction(Vec positions, PetscScalar *objectiveFunction)
{
    PetscScalar unused;
    PetscCall(
        physics->ComputeObjectiveConstraints(objectiveFunction,
                                            &unused,
                                            positions,
                                            opt->Emin,
                                            opt->Emax,
                                            opt->penal,
                                            opt->volfrac);
    );

    return 0;
}
