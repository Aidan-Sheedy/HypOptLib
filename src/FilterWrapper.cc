#include "FilterWrapper.h"

PetscErrorCode FilterWrapper::filterDesignVariable(Vec unfiltered, Vec filtered)
{
    Vec unused;
    PetscCall(VecDuplicate(unfiltered, &unused));

    this->filter->FilterProject(unfiltered, filtered, unused, PETSC_FALSE, 0, 0);

    VecDestroy(&unused);

    return 0;
}


PetscErrorCode FilterWrapper::filterSensitivities(Vec position, Vec sensitivities, Vec* constraintSensitivities)
{
    // Used: x, dfdx, dgdx
    //       position, sensitivities, constraintSensitivities
    Vec unused;
    PetscCall(VecDuplicate(position, &unused));


    PetscCall(
        filter->Gradients(position, unused, sensitivities, 0, constraintSensitivities, PETSC_FALSE, 0, 0)
        );
    
    VecDestroy(&unused);

    return 0;
}
