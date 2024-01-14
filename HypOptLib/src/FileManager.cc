/***************************************************************************//**
 * @file FileManager.cc
 *
 * This file implements the FileManager class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#include "FileManager.h"

#include "PetscExtensions.h"
#include "HypOptException.h"
#include <petscviewerhdf5.h>

#include <limits.h>

/**
 * File-related constant attributes
 *
 * This can be defined in the header in C++17 using the inline keyword,
 * but there is no guarantee that all the user's compiler suports C++17.
 **/

/* HDF5 groups */
const std::string FileManager::stateGroup    = "/Dataset/State";
const std::string FileManager::dataGroup     = "/Dataset";
const std::string FileManager::settingsGroup = "/Setting";

/* System parameters and settings */
const std::string FileManager::volfracName   = "volfrac";
const std::string FileManager::timestepName  = "dt";
const std::string FileManager::tempName      = "T";
const std::string FileManager::xDimName      = "nelx";
const std::string FileManager::yDimName      = "nely";
const std::string FileManager::zDimName      = "nelz";
const std::string FileManager::penaltyName   = "penal";
const std::string FileManager::minRadiusName = "rmin";
const std::string FileManager::numStepName   = "stepN";
const std::string FileManager::numSampleName = "sampleN";
const std::string FileManager::chainOrdName  = "ch_order";

/* Simulation Data */
const std::string FileManager::hamiltonianName           = "Hamiltonian";
const std::string FileManager::complainceName            = "Compliance";
const std::string FileManager::temperatureName           = "Temperature";
const std::string FileManager::lagrangianMultipliarName  = "Lambda";
const std::string FileManager::iterationComputeTimeName  = "Iteration Compute Time";
const std::string FileManager::timestepsName             = "Timestep";
const std::string FileManager::energyErrorsName          = "Energy Error";
const std::string FileManager::finalPositionName         = "Final Position";
const std::string FileManager::finalVelocityName         = "Final Velocity";
const std::string FileManager::finalEvenNHPositionName   = "Final even NH Pos";
const std::string FileManager::finalEvenNHVelocityName   = "Final even NH Vel";
const std::string FileManager::finalOddNHPositionName    = "Final odd NH Pos";
const std::string FileManager::finalOddNHVelocityName    = "Final odd NH Vel";
const std::string FileManager::finalStateFieldName       = "Final State Field";

/* Test file names */
const std::string FileManager::initialPositionName = "Initial Position";
const std::string FileManager::initialVelocityName = "Initial Velocity";

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
                                           std::string  filePath,
                                           bool         randomStartingValues,
                                           std::string  initialConditionsFile)
{
    PetscErrorCode errorStatus = 0;

    if (".h5" != filePath.substr(filePath.size() - 3))
    {
        throw HypOptException("Bad input file path. Must have file extension \'.h5\'!");
    }

    /* Avoid overwriting existing files */
    this->saveFilePath = autoAppendFilePath(filePath.substr(0, filePath.size() - 3), ".h5");

    PetscPrintf(PETSC_COMM_WORLD, "# Saving to file: %s\n", saveFilePath.c_str());

    /* Open and set up the save file */
    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, saveFilePath.c_str(), FILE_MODE_WRITE, &(saveFileHDF5)));
    PetscCall(PetscViewerHDF5WriteGroup(saveFileHDF5, settingsGroup.c_str()));
    PetscCall(PetscViewerHDF5WriteGroup(saveFileHDF5, this->stateGroup.c_str()));

    /* Write all settings */
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), volfracName.c_str(),      PETSC_SCALAR,   &volfrac));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), timestepName.c_str(),     PETSC_SCALAR,   &timestep));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), tempName.c_str(),         PETSC_SCALAR,   &temperature));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), xDimName.c_str(),         PETSC_INT,      &numberElementsX));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), yDimName.c_str(),         PETSC_INT,      &numberElementsY));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), zDimName.c_str(),         PETSC_INT,      &numberElementsZ));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), penaltyName.c_str(),      PETSC_SCALAR,   &penalty));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), minRadiusName.c_str(),    PETSC_SCALAR,   &filterRadius));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), numStepName.c_str(),      PETSC_INT,      &numberSteps));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), numSampleName.c_str(),    PETSC_INT,      &numberSamples)); /** @todo IMPLEMENT */
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), chainOrdName.c_str(),     PETSC_INT,      &NoseHooverChainOrder));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), "Random Init Conditions", PETSC_BOOL,     &randomStartingValues));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), "Init Conditions File",   PETSC_STRING,   initialConditionsFile.c_str()));

    PetscCall(PetscViewerDestroy(&(saveFileHDF5)));

    return errorStatus;
}

PetscErrorCode FileManager::saveIteration(PetscInt iteration, Vec positions)
{
    PetscErrorCode errorStatus = 0;

    /* Open the file */
    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, this->saveFilePath.c_str(), FILE_MODE_APPEND, &(saveFileHDF5)));

    /* Set positions name */
    std::string stateName = "iteration";
    stateName.append(std::to_string(iteration));

    /* Save the iteration */
    PetscCall(PetscObjectSetName((PetscObject)positions, stateName.c_str()));
    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, stateGroup.c_str()));
    PetscCall(VecView(positions, saveFileHDF5));

    /* Cleanup */
    PetscCall(PetscViewerHDF5PopGroup(saveFileHDF5));
    PetscCall(PetscViewerDestroy(&saveFileHDF5));

    return errorStatus;
}

PetscErrorCode FileManager::saveFinalState( bool                     saveHamiltionan,
                                            HypOptParameters         finalState,
                                            std::vector<PetscScalar> hamiltonians,
                                            std::vector<PetscScalar> compliance,
                                            std::vector<PetscScalar> temperatures,
                                            std::vector<PetscScalar> LagrangeMultipliers,
                                            std::vector<PetscScalar> iterationTimes,
                                            std::vector<PetscScalar> timesteps,
                                            std::vector<PetscScalar> energyErrors,
                                            std::vector<PetscScalar> volFracs)
{
    PetscErrorCode errorStatus = 0;

    /* Set up save file */
    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, saveFilePath.c_str(), FILE_MODE_APPEND, &(saveFileHDF5)));
    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, dataGroup.c_str()));

    PetscPrintf(PETSC_COMM_WORLD, "# Saving final state...");

    /* Save std vectors */
    if (saveHamiltionan)
    {
        PetscCall(HDF5SaveStdVector(saveFileHDF5, hamiltonians, hamiltonianName.c_str()));
        PetscCall(HDF5SaveStdVector(saveFileHDF5, compliance,   complainceName.c_str()));
    }
    PetscCall(HDF5SaveStdVector(saveFileHDF5, temperatures,          temperatureName.c_str()));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, LagrangeMultipliers,   lagrangianMultipliarName.c_str()));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, iterationTimes,        iterationComputeTimeName.c_str()));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, timesteps,             timestepsName.c_str()));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, energyErrors,          energyErrorsName.c_str()));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, volFracs,          "Volume Fraction"));

    /* save Petsc type vectors */
    PetscCall(PetscObjectSetName((PetscObject)(finalState.position),                finalPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.velocity),                finalVelocityName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.evenNoseHooverPosition),  finalEvenNHPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.evenNoseHooverVelocity),  finalEvenNHVelocityName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.oddNoseHooverPosition),   finalOddNHPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)(finalState.oddNoseHooverVelocity),   finalOddNHVelocityName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)(physics->GetStateField()),           finalStateFieldName.c_str()));

    PetscCall(VecView(finalState.position,                  saveFileHDF5));
    PetscCall(VecView(finalState.velocity,                  saveFileHDF5));
    PetscCall(VecView(finalState.evenNoseHooverPosition,    saveFileHDF5));
    PetscCall(VecView(finalState.evenNoseHooverVelocity,    saveFileHDF5));
    PetscCall(VecView(finalState.oddNoseHooverPosition,     saveFileHDF5));
    PetscCall(VecView(finalState.oddNoseHooverVelocity,     saveFileHDF5));
    PetscCall(VecView(physics->GetStateField(),             saveFileHDF5));

    PetscPrintf(PETSC_COMM_WORLD, "# ...saved!\n");

    /* Close file */
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

PetscErrorCode FileManager::getFinalStateVectors(   std::string         filePath,
                                                    HypOptParameters    finalState,
                                                    Vec                 finalStateField)
{
    if (!doesFileExist(filePath))
    {
        throw HypOptException("Invalid restart file, please check your file path to make sure it exists.");
    }

    /* Can't assume that the provided vectors are already correctly named, so rename them. */
    PetscCall(PetscObjectSetName((PetscObject)finalState.position,               finalPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)finalState.velocity,               finalVelocityName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)finalState.evenNoseHooverPosition, finalEvenNHPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)finalState.evenNoseHooverVelocity, finalEvenNHVelocityName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)finalState.oddNoseHooverPosition,  finalOddNHPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)finalState.oddNoseHooverVelocity,  finalOddNHVelocityName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)finalStateField,                   finalStateFieldName.c_str()));

    /* Get the vectors from the appropriate file. */
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalState.position));
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalState.velocity));
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalState.evenNoseHooverPosition));
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalState.evenNoseHooverVelocity));
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalState.oddNoseHooverPosition));
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalState.oddNoseHooverVelocity));
    PetscCall(HDF5GetSavedVec(filePath, dataGroup, finalStateField));

    return 0;
}

PetscErrorCode FileManager::HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<PetscScalar> vector, const char * vectorName)
{
    PetscErrorCode errorStatus = 0;

    /* Create vectors */
    Vec outputVector;
    PetscCall(VecCreate(PETSC_COMM_WORLD, &outputVector));
    PetscCall(VecSetSizes(outputVector, PETSC_DECIDE, vector.size()));
    PetscCall(VecSetFromOptions(outputVector));

    /* Get values from standard vectors */
    PetscCall(PetscExtensions::VecParallelFromStdVector(vector, outputVector));
    PetscCall(PetscObjectSetName((PetscObject)outputVector, vectorName));

    /* Save Vectors */
    PetscCall(VecView(outputVector, HDF5saveFile));
    PetscCall(VecDestroy(&outputVector));

    return errorStatus;
}

PetscErrorCode FileManager::HDF5GetSavedVec(std::string filePath, std::string location, Vec vector)
{
    PetscErrorCode errorStatus = 0;

    /* Open the file */
    PetscViewer saveFile;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &saveFile));

    /* Grab the vector*/
    PetscCall(PetscViewerHDF5PushGroup(saveFile, location.c_str()));
    PetscCall(VecLoad(vector, saveFile));

    /* Cleanup */
    PetscCall(PetscViewerHDF5PopGroup(saveFile));
    PetscCall(PetscViewerDestroy(&saveFile));

    return errorStatus;
}

PetscErrorCode FileManager::getHDF5Settings(std::string             filePath,
                                            PetscInt               *noseHooverChainOrder,
                                            PetscScalar            *volumeFraction,
                                            PetscScalar            *timestep,
                                            PetscScalar            *targetTemperature,
                                            PetscScalar            *penalty,
                                            PetscScalar            *minimumFilterRadius,
                                            std::vector<uint32_t>  *gridDimensions)
{
    PetscInt    numberElementsX;
    PetscInt    numberElementsY;
    PetscInt    numberElementsZ;
    PetscViewer saveFile;

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &saveFile));

    /* Retrieve all parameters */
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), chainOrdName.c_str(),   PETSC_INT,      NULL, noseHooverChainOrder));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), volfracName.c_str(),    PETSC_SCALAR,   NULL, volumeFraction));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), timestepName.c_str(),   PETSC_SCALAR,   NULL, timestep));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), tempName.c_str(),       PETSC_SCALAR,   NULL, targetTemperature));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), penaltyName.c_str(),    PETSC_SCALAR,   NULL, penalty));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), minRadiusName.c_str(),  PETSC_SCALAR,   NULL, minimumFilterRadius));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), xDimName.c_str(),       PETSC_INT,      NULL, &numberElementsX));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), yDimName.c_str(),       PETSC_INT,      NULL, &numberElementsY));
    PetscCall(PetscViewerHDF5ReadAttribute(saveFile, settingsGroup.c_str(), zDimName.c_str(),       PETSC_INT,      NULL, &numberElementsZ));

    /* Convert to vector */
    gridDimensions->push_back(numberElementsX);
    gridDimensions->push_back(numberElementsY);
    gridDimensions->push_back(numberElementsZ);

    /* Cleanup */
    PetscCall(PetscViewerDestroy(&(saveFile)));

    return 0;
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
            if (INT_MAX == fileNumber)
            {
                throw HypOptException("Can not find a unique file for the name provided. You have.... so many files...");
            }
            suffix = " (" + std::to_string(++fileNumber) + ")";
        }
    } while (doesFileExist(filePath));

    return filePath;
}

PetscErrorCode FileManager::saveInitialConditions(Vec positions, Vec velocities, std::string filePath)
{
    if (".h5" != filePath.substr(filePath.size() - 3))
    {
        throw HypOptException("Bad input file path. Must have file extension \'.h5\'!");
    }

    /* Open and set up the save file */
    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &(saveFileHDF5)));

    /* Can't assume that the provided vectors are already correctly named, so rename them. */
    PetscCall(PetscObjectSetName((PetscObject)positions,  initialPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)velocities, initialVelocityName.c_str()));

    /* Get the vectors from the appropriate file. */
    PetscCall(VecView(positions,  saveFileHDF5));
    PetscCall(VecView(velocities, saveFileHDF5));

    PetscCall(PetscViewerDestroy(&saveFileHDF5));

    return 0;
}

PetscErrorCode FileManager::loadInitialConditions(Vec positions, Vec velocities, std::string filePath)
{
    if (".h5" != filePath.substr(filePath.size() - 3))
    {
        throw HypOptException("Bad input file path. Must have file extension \'.h5\'!");
    }

    /* Can't assume that the provided vectors are already correctly named, so rename them. */
    PetscCall(PetscObjectSetName((PetscObject)positions,  initialPositionName.c_str()));
    PetscCall(PetscObjectSetName((PetscObject)velocities, initialVelocityName.c_str()));

    PetscCall(HDF5GetSavedVec(filePath, "", positions));
    PetscCall(HDF5GetSavedVec(filePath, "", velocities));

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