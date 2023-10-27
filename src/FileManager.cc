/**
 * @file Hyperoptimization.h
 * 
 * This file describes the layout of the Hyperoptimization class.
 * 
 * @todo LICENSE
**/

#include "FileManager.h"

#include "PetscExtensions.h"
#include "HypOptException.h"
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

PetscErrorCode FileManager::HDF5GetSavedVec(std::string filePath, std::string location, Vec vector)
{
    PetscErrorCode errorStatus = 0;

    PetscViewer saveFile;


    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &saveFile));

    PetscCall(PetscViewerHDF5PushGroup(saveFile, location.c_str()));
    PetscCall(VecLoad(vector, saveFile));
    PetscCall(PetscViewerHDF5PopGroup(saveFile));

    PetscCall(PetscViewerDestroy(&saveFile));

    return errorStatus;
}

bool FileManager::doesFileExist(std::string filePath)
{
    std::ifstream testFile;
    testFile.open(filePath);
    if (!testFile)
    {
        return false;
    }
    testFile.close();

    return true;
}

std::string FileManager::autoAppendFilePath(std::string fileName, std::string fileExtension)
{
    std::string filePath;
    std::string suffix = "";
    PetscInt    fileNumber = 0;

    do
    {
        filePath = fileName + suffix + fileExtension;
        if (doesFileExist(filePath))
        {
            suffix = " (" + std::to_string(++fileNumber) + ")";
        }
    } while (doesFileExist(filePath));

    return filePath;
}

PetscErrorCode FileManager::initializeHDF5(PetscScalar  volfrac,
                                           PetscScalar  timestep,
                                           PetscScalar  temperature,
                                           PetscInt     numberElementsX,
                                           PetscInt     numberElementsY,
                                           PetscInt     numberElementsZ,
                                           PetscScalar  penalty,
                                           PetscScalar  filterRadius,
                                           PetscInt     numberSteps,
                                           PetscInt     numberSamples,
                                           PetscInt     NoseHooverChainOrder,
                                           std::string  filePath)
{
    PetscErrorCode errorStatus = 0;

    if ("" == filePath)
    {
        filePath = autoAppendFilePath("hypopt_output", ".h5");
    }
    else
    {
        if (".h5" != filePath.substr(filePath.size() - 3))
        {
            throw HypOptException("Bad input file path. Must have file extension \'.h5!");
        }
        filePath = autoAppendFilePath(filePath.substr(0, filePath.size() - 3), ".h5");
    }

    this->saveFilePath = filePath;

    PetscPrintf(PETSC_COMM_WORLD, "# Saving to file: %s\n", filePath.c_str());

    PetscViewer saveFileHDF5;

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &(saveFileHDF5)));

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
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "ch_order",   PETSC_INT,      &NoseHooverChainOrder));

    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, dataGroup.c_str()));

    /* Save Mesh Settings for XDMF later */
    Vec xDim;
    Vec yDim;
    Vec zDim;

    generateDimensionVector(&xDim, numberElementsX+1);
    generateDimensionVector(&yDim, numberElementsY+1);
    generateDimensionVector(&zDim, numberElementsZ+1);

    PetscCall(PetscObjectSetName((PetscObject)(xDim), "xMesh"));
    PetscCall(PetscObjectSetName((PetscObject)(yDim), "yMesh"));
    PetscCall(PetscObjectSetName((PetscObject)(zDim), "zMesh"));

    PetscCall(VecView(xDim, saveFileHDF5));
    PetscCall(VecView(yDim, saveFileHDF5));
    PetscCall(VecView(zDim, saveFileHDF5));



    PetscCall(PetscViewerDestroy(&(saveFileHDF5)));
    PetscCall(VecDestroy(&xDim));
    PetscCall(VecDestroy(&yDim));
    PetscCall(VecDestroy(&zDim));

    return errorStatus;
}

PetscErrorCode FileManager::generateDimensionVector(Vec *dimension, PetscInt meshDimension)
{
    PetscCall(VecCreate(PETSC_COMM_WORLD, dimension));

    PetscCall(VecSetSizes(*dimension, PETSC_DETERMINE, meshDimension));

    PetscCall(VecSetFromOptions(*dimension));

    PetscScalar *dimVec;

    PetscCall(VecGetArray(*dimension, &dimVec));

    PetscInt xLowOwnership;
    PetscInt xHighOwnership;

    PetscCall(VecGetOwnershipRange(*dimension, &xLowOwnership, &xHighOwnership));

    for (PetscInt i = 0; i < (xHighOwnership - xLowOwnership); i++)
    {
        dimVec[i] = (i+(PetscScalar)xLowOwnership)/(PetscScalar)meshDimension;
    }

    PetscCall(VecRestoreArray(*dimension, &dimVec));

    return 0;
}

// PetscErrorCode FileManager::preparetoIterate(Vec *positions, PetscInt firstTimestep)
// {
//     PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, this->saveFilePath.c_str(), FILE_MODE_APPEND, &(iterationViewer)));
//     PetscCall(PetscViewerHDF5PushTimestepping(iterationViewer));
//     PetscCall(PetscViewerHDF5PushGroup(iterationViewer, stateGroup.c_str()));
//     PetscCall(PetscObjectSetName((PetscObject)*positions, "Positions"));

//     /* If VecCreate is used to make generate the iterating vector, then this is not necessary. */
//     PetscCall(DMSetOutputSequenceNumber(physics->GetDM(), 0, firstTimestep));
//     return 0;
// }

// PetscErrorCode FileManager::finishIterating()
// {
//     PetscCall(PetscViewerHDF5PopGroup(iterationViewer));
//     PetscCall(PetscViewerHDF5PopTimestepping(iterationViewer));
//     PetscCall(PetscViewerDestroy(&iterationViewer));
//     return 0;
// }

// PetscErrorCode FileManager::saveIteration(PetscInt iteration, Vec positions, PetscInt nextTimestep)
// {
//     PetscErrorCode errorStatus = 0;
//     PetscCall(PetscViewerHDF5SetTimestep(iterationViewer, nextTimestep));
//     PetscCall(VecView(positions, iterationViewer));

//      /* If VecCreate is used to make generate the iterating vector, then this use
//       * PetscViewerHDF5IncrementTimestep instead. 
//       */
//     PetscCall(DMSetOutputSequenceNumber(physics->GetDM(), iteration+1, nextTimestep));
//     return errorStatus;
// }

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

PetscErrorCode FileManager::getHDF5Settings(std::string             filePath,
                                            PetscInt               *noseHooverChainOrder,
                                            PetscInt               *numSavedIterations,
                                            PetscScalar            *volumeFraction,
                                            PetscScalar            *timestep,
                                            PetscScalar            *targetTemperature,
                                            PetscScalar            *penalty,
                                            PetscScalar            *minimumFilterRadius,
                                            std::vector<uint32_t>  *gridDimensions)
{
    /* Open file */
    PetscViewer saveFile;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &saveFile));

    /* Start with simple values. No default values are used, if the value is not found the file is invalid. */

    PetscInt numberElementsX;
    PetscInt numberElementsY;
    PetscInt numberElementsZ;

    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "ch_order",    PETSC_INT,      NULL, noseHooverChainOrder));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "stepN",       PETSC_INT,      NULL, numSavedIterations));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "volfrac",     PETSC_SCALAR,   NULL, volumeFraction));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "dt",          PETSC_SCALAR,   NULL, timestep));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "T",           PETSC_SCALAR,   NULL, targetTemperature));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "penal",       PETSC_SCALAR,   NULL, penalty));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "rmin",        PETSC_SCALAR,   NULL, minimumFilterRadius));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "nelx",        PETSC_INT,      NULL, &numberElementsX));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "nely",        PETSC_INT,      NULL, &numberElementsY));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, "/Setting", "nelz",        PETSC_INT,      NULL, &numberElementsZ));

    gridDimensions->push_back(numberElementsX);
    gridDimensions->push_back(numberElementsY);
    gridDimensions->push_back(numberElementsZ);

    PetscCall(PetscViewerDestroy(&(saveFile)));

    return 0;
}

PetscErrorCode FileManager::getHDF5Vectors( std::string         filePath,
                                            HypOptParameters    finalState,
                                            Vec                 finalStateField)
{
    /* Use number of saved iterations to get final position */
    PetscCall(PetscObjectSetName((PetscObject)finalState.position,               "Final Position"));
    PetscCall(PetscObjectSetName((PetscObject)finalState.velocity,               "Final Velocity"));
    PetscCall(PetscObjectSetName((PetscObject)finalState.evenNoseHooverPosition, "Final even NH Pos"));
    PetscCall(PetscObjectSetName((PetscObject)finalState.evenNoseHooverVelocity, "Final even NH Vel"));
    PetscCall(PetscObjectSetName((PetscObject)finalState.oddNoseHooverPosition,  "Final odd NH Pos"));
    PetscCall(PetscObjectSetName((PetscObject)finalState.oddNoseHooverVelocity,  "Final odd NH Vel"));
    PetscCall(PetscObjectSetName((PetscObject)finalStateField,                   "Final State Field"));

    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalState.position));
    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalState.velocity));
    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalState.evenNoseHooverPosition));
    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalState.evenNoseHooverVelocity));
    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalState.oddNoseHooverPosition));
    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalState.oddNoseHooverVelocity));
    PetscCall(HDF5GetSavedVec(filePath, "/Dataset", finalStateField));

    return 0;
}

PetscErrorCode FileManager::saveFinalState( bool saveHamiltionan,
                                            HypOptParameters finalState,
                                            std::vector<PetscScalar> hamiltonians,
                                            std::vector<PetscScalar> compliance,
                                            std::vector<PetscScalar> temperatures,
                                            std::vector<PetscScalar> LagrangeMultipliers,
                                            std::vector<PetscScalar> iterationTimes)
{
    PetscErrorCode errorStatus = 0;

    /* Set up save file */
    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, saveFilePath.c_str(), FILE_MODE_APPEND, &(saveFileHDF5)));
    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, dataGroup.c_str()));

    /* Save std vectors */
    if (saveHamiltionan)
    {
        PetscCall(HDF5SaveStdVector(saveFileHDF5, hamiltonians, "Hamiltonian"));
        PetscCall(HDF5SaveStdVector(saveFileHDF5, compliance,   "Compliance"));
    }
    PetscCall(HDF5SaveStdVector(saveFileHDF5, temperatures,          "Temperature"));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, LagrangeMultipliers,   "Lambda"));
    // PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->genericData,           "Volume Fraction"));
    // PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->genericData2,          "Max Position"));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, iterationTimes,        "Iteration Compute Time"));

    PetscPrintf(PETSC_COMM_WORLD, "# Saving final state...");

    /* save Petsc type vectors */
    PetscCall(PetscObjectSetName((PetscObject)(finalState.position),                "Final Position"));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.velocity),                "Final Velocity"));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.evenNoseHooverPosition),  "Final even NH Pos"));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.evenNoseHooverVelocity),  "Final even NH Vel"));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.oddNoseHooverPosition),   "Final odd NH Pos"));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.oddNoseHooverVelocity),   "Final odd NH Vel"));
    PetscCall(PetscObjectSetName((PetscObject)(physics->GetStateField()),           "Final State Field"));

    PetscCall(VecView(finalState.position,                  saveFileHDF5));
    PetscCall(VecView(finalState.velocity,                  saveFileHDF5));
    PetscCall(VecView(finalState.evenNoseHooverPosition,    saveFileHDF5));
    PetscCall(VecView(finalState.evenNoseHooverVelocity,    saveFileHDF5));
    PetscCall(VecView(finalState.oddNoseHooverPosition,     saveFileHDF5));
    PetscCall(VecView(finalState.oddNoseHooverVelocity,     saveFileHDF5));
    PetscCall(VecView(physics->GetStateField(),             saveFileHDF5));

    PetscPrintf(PETSC_COMM_WORLD, "...saved!\n");

    /* Close file */
    PetscCall(PetscViewerHDF5PopGroup(saveFileHDF5));
    PetscCall(PetscViewerDestroy(&saveFileHDF5));

    return errorStatus;
}
