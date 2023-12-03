#include "SensitivitiesWrapper.h"

#include "PetscExtensions.h"

/** @todo figure out the opt-> issue (probably just pass these in?)*/
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
                                      volfrac);
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
                                            volfrac);
    );

    return 0;
}
