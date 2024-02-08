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
                                      Emin,
                                      Emax,
                                      penal,
                                      volfrac,
                                      &petscSolverIterationCount);
    );

    return 0;
}

PetscErrorCode SensitivitiesWrapper::computeObjectiveFunction(Vec filteredPositions, PetscScalar *objectiveFunction)
{
    PetscScalar unused;
    PetscCall(
        physics->ComputeObjectiveConstraints(objectiveFunction,
                                            &unused,
                                            filteredPositions,
                                            Emin,
                                            Emax,
                                            penal,
                                            volfrac,
                                            &petscSolverIterationCount);
    );

    return 0;
}
