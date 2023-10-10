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

PetscErrorCode FileManager::initializeHDF5(PetscScalar volfrac,
                                           PetscScalar timestep,
                                           PetscScalar temperature,
                                           PetscInt numberElementsX,
                                           PetscInt numberElementsY,
                                           PetscInt numberElementsZ,
                                           PetscScalar penalty,
                                           PetscScalar filterRadius,
                                           PetscInt numberSteps,
                                           PetscInt numberSamples,
                                           PetscScalar NoseHooverChainOrder)
{
    PetscPrintf(PETSC_COMM_WORLD, "\t- Entered File Manager\n");
    PetscErrorCode errorStatus = 0;

    std::string filename = "hypopt_output";
    std::string fileExtension = ".h5";
    std::string fullPath = "";

    // Check PETSc input for a work directory
    // char      filenameChar[PETSC_MAX_PATH_LEN];
    // PetscBool flg = PETSC_FALSE;
    // PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);

    // If input, change path of the file in filename
    // if (flg) {
    //     fullPath.append(filenameChar);
    //     fullPath.append("/");
    // }
    // else
    // {
    //     fullPath.append("output/");
    // }
    fullPath.append(filename);

    // std::filesystem::path filePath{fullPath + fileExtension};
    // int i = 0;
    std::string appendor = "";
    // while (std::filesystem::exists(filePath))
    // {
    //     i++;
    //     appendor = std::string(i);
    //     filePath = std::filesystem::path(fullPath + appendor + fileExtension);
    // }

    fullPath.append(appendor + fileExtension);

    this->saveFilePath = fullPath;

    PetscPrintf(PETSC_COMM_WORLD, "Trying to save to file: %s\n", fullPath.c_str());

    PetscViewer saveFileHDF5;

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, fullPath.c_str(), FILE_MODE_WRITE, &(saveFileHDF5)));

    // std::time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    // std::string timeString = std::asctime(currentTime);
    // PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/", "create_date", PETSC_STRING, &()));

    PetscCall(PetscViewerHDF5WriteGroup(saveFileHDF5, "/Setting"));
    PetscCall(PetscViewerHDF5WriteGroup(saveFileHDF5, this->stateGroup.c_str()));

    /* Write all settings */
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "volfrac",    PETSC_SCALAR,   &volfrac));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "dt",         PETSC_SCALAR,   &timestep));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "T",          PETSC_SCALAR,   &temperature));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "nelx",       PETSC_INT,      &numberElementsX));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "nely",       PETSC_INT,      &numberElementsY));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "nelz",       PETSC_INT,      &numberElementsZ));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "penal",      PETSC_SCALAR,   &penalty));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "rmin",       PETSC_SCALAR,   &filterRadius));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "stepN",      PETSC_INT,      &numberSteps));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "sampleN",    PETSC_INT,      &numberSamples)); /** @todo IMPLEMENT */
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "ch_Order",   PETSC_SCALAR,   &NoseHooverChainOrder));

    PetscCall(PetscViewerDestroy(&(saveFileHDF5)));

    return errorStatus;
}

PetscErrorCode FileManager::saveIteration(PetscInt iteration, Vec positions)
{
    PetscErrorCode errorStatus = 0;

    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, this->saveFilePath.c_str(), FILE_MODE_APPEND, &(saveFileHDF5)));

    // Set positions name
    std::string stateName = "iteration";
    stateName.append(std::to_string(iteration));

    PetscCall(PetscObjectSetName((PetscObject)positions, stateName.c_str()));
    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, stateGroup.c_str()));
    PetscCall(VecView(positions, saveFileHDF5));
    PetscCall(PetscViewerHDF5PopGroup(saveFileHDF5));

    PetscCall(PetscViewerDestroy(&saveFileHDF5));

    return errorStatus;
}

std::string FileManager::getSaveFilePath()
{
    return saveFilePath;
}

std::string FileManager::getDataGroup()
{
    return dataGroup;
}
