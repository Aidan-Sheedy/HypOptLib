/**
 * @file Hyperoptimization.h
 * 
 * This file describes the layout of the Hyperoptimization class.
 * 
 * @todo LICENSE
**/

#include "FileManager.h"

#include "PetscExtensions.h"
#include <petscviewerhdf5.h>

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

PetscErrorCode FileManager::HDF5GetSavedVec(std::string filePath, std::string location, Vec *vector)
{
    PetscErrorCode errorStatus = 0;

    PetscViewer saveFile;

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &saveFile));

    PetscCall(PetscViewerHDF5PushGroup(saveFile, location.c_str()));
    PetscCall(VecLoad(*vector, saveFile));
    PetscCall(PetscViewerHDF5PopGroup(saveFile));

    PetscCall(PetscViewerDestroy(&saveFile));

    return errorStatus;
}