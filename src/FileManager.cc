/**
 * @file Hyperoptimization.h
 * 
 * This file describes the layout of the Hyperoptimization class.
 * 
 * @todo LICENSE
**/

#include "FileManager.h"

#include "PetscExtensions.h"

PetscErrorCode FileManager::HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<PetscScalar> vector, const char * vectorName)
{
    PetscErrorCode errorStatus = 0;

    /* Create vectors */
    Vec outputVector;

    PetscCall(VecCreate(PETSC_COMM_WORLD, &outputVector));
    PetscCall(VecSetSizes(outputVector, PETSC_DECIDE, vector.size()));
    PetscCall(VecSetFromOptions(outputVector));

    /* Get values from standard vectors */
    PetscCall(PetscExtensions::VecParallelFromStdVector(outputVector, vector));
    PetscCall(PetscObjectSetName((PetscObject)outputVector, vectorName));

    /* Save Vectors */
    PetscCall(VecView(outputVector, HDF5saveFile));
    PetscCall(VecDestroy(&outputVector));

    return errorStatus;
}