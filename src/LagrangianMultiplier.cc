
#include "LagrangianMultiplier.h"

PetscScalar LagrangianMultiplier::computeLagrangianMultiplier(Vec positions, Vec C, PetscInt numParticles, PetscScalar *returnValue)
{
    PetscErrorCode errorStatus = 0;

    // Copies of input parameters
    Vec positionCopy;
    Vec CCopy;

    PetscCall(VecDuplicate(positions, &positionCopy));
    PetscCall(VecDuplicate(C, &CCopy));
    PetscCall(VecCopy(positions, positionCopy));
    PetscCall(VecCopy(C, CCopy));

    // Filter positions
    PetscCall(VecCopy(positionCopy, opt->x));
    PetscCall(this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta));
    PetscCall(VecCopy(opt->xTilde, positionCopy));

    // Filter positions
    PetscCall(VecCopy(CCopy, opt->x));
    PetscCall(this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta));
    PetscCall(VecCopy(opt->xTilde, CCopy));

    // Sum values
    PetscScalar positionSum;
    PetscCall(VecSum(positionCopy, &positionSum));
    PetscScalar CSum;
    PetscCall(VecSum(CCopy, &CSum));


    *returnValue = (numParticles * opt->volfrac - positionSum )/ CSum;

    PetscCall(VecDestroy(&positionCopy));
    PetscCall(VecDestroy(&CCopy));

    return errorStatus;
}