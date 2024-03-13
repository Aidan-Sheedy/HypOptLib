/***************************************************************************//**
 * @file FileManager.cc
 *
 * This file implements the FileManager class.
 *
 * @author Aidan Sheedy
 *
 * Copyright (C) 2024 Aidan Sheedy
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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
const std::string FileManager::saveFreqName = "saveFreq";
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
                                           PetscInt     saveFrequency,
                                           PetscInt     NoseHooverChainOrder,
                                           std::string  filePath,
                                           bool         randomStartingValues,
                                           std::string  initialConditionsFile,
                                           DomainCoordinates domain,
                                           std::vector<BoundaryCondition> boundaryConditions)
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
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), saveFreqName.c_str(),     PETSC_INT,      &saveFrequency));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), chainOrdName.c_str(),     PETSC_INT,      &NoseHooverChainOrder));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), "Random Init Conditions", PETSC_BOOL,     &randomStartingValues));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, settingsGroup.c_str(), "Init Conditions File",   PETSC_STRING,   initialConditionsFile.c_str()));

    /* Mesh and Boundary Conditions */
    saveBoundaryConditions(saveFileHDF5, boundaryConditions);
    saveDomainInformation(saveFileHDF5, domain);

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
                                            std::vector<PetscScalar> volFracs,
                                            std::vector<PetscInt>    solverIterationsSensitivity,
                                            std::vector<PetscInt>    solverIterationsHamiltonian)
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
    PetscCall(HDF5SaveStdVector(saveFileHDF5, solverIterationsSensitivity, "FEA Solver Itr Sensitivity"));
    PetscCall(HDF5SaveStdVector(saveFileHDF5, solverIterationsHamiltonian, "FEA Solver Itr Hamiltonian"));

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

    PetscPrintf(PETSC_COMM_WORLD, "...saved!\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Output to file: %s\n", saveFilePath.c_str());

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

template<typename T> PetscErrorCode FileManager::HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<T> vector, const char * vectorName)
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
                                            std::vector<uint32_t>  *gridDimensions,
                                            DomainCoordinates      *domain,
                                            std::vector<BoundaryCondition> *boundaryConditions)
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

    /* Get mesh and boundary condition info */
    getBoundaryConditions(saveFile, boundaryConditions);
    getDomainInformation(saveFile, domain);

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

PetscErrorCode FileManager::saveBoundaryConditions(PetscViewer saveFileHDF5, std::vector<BoundaryCondition> boundaryConditions)
{
    PetscBool saveTrue = PETSC_TRUE;
    PetscBool saveFalse = PETSC_FALSE;

    uint32_t boundaryCondtionNumber = 0;
    for (auto boundaryCondition : boundaryConditions)
    {
        std::string parent = settingsGroup + "/boundaryConditions/boundaryCondition_" + std::to_string(boundaryCondtionNumber);
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "type",       PETSC_INT,      &boundaryCondition.type));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "xRangeMin",  PETSC_SCALAR,   &boundaryCondition.xRange[0]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "xRangeMax",  PETSC_SCALAR,   &boundaryCondition.xRange[1]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "yRangeMin",  PETSC_SCALAR,   &boundaryCondition.yRange[0]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "yRangeMax",  PETSC_SCALAR,   &boundaryCondition.yRange[1]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "zRangeMin",  PETSC_SCALAR,   &boundaryCondition.zRange[0]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "zRangeMax",  PETSC_SCALAR,   &boundaryCondition.zRange[1]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "zRangeMax",  PETSC_SCALAR,   &boundaryCondition.zRange[1]));
        PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "value",      PETSC_SCALAR,   &boundaryCondition.value));

        if (boundaryCondition.degreesOfFreedom.find(0) != boundaryCondition.degreesOfFreedom.end())
            PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "dofX",  PETSC_BOOL, &saveTrue));
        else
            PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "dofX",  PETSC_BOOL, &saveFalse));

        if (boundaryCondition.degreesOfFreedom.find(1) != boundaryCondition.degreesOfFreedom.end())
            PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "dofY",  PETSC_BOOL, &saveTrue));
        else
            PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "dofY",  PETSC_BOOL, &saveFalse));

        if (boundaryCondition.degreesOfFreedom.find(2) != boundaryCondition.degreesOfFreedom.end())
            PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "dofZ",  PETSC_BOOL, &saveTrue));
        else
            PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "dofZ",  PETSC_BOOL, &saveFalse));

        boundaryCondtionNumber++;
    }

    return 0;
}

PetscErrorCode FileManager::getBoundaryConditions(PetscViewer saveFileHDF5, std::vector<BoundaryCondition> *boundaryConditions)
{
    PetscBool dofX;
    PetscBool dofY;
    PetscBool dofZ;
    uint32_t boundaryCondtionNumber = 0;
    PetscBool hasDataset = PETSC_FALSE;
    std::string parent = settingsGroup + "/boundaryConditions/boundaryCondition_" + std::to_string(boundaryCondtionNumber);

    PetscCall(PetscViewerHDF5HasGroup(saveFileHDF5, parent.c_str(), &hasDataset));
    while (hasDataset)
    {
        boundaryConditions->push_back(BoundaryCondition());
        boundaryConditions->at(boundaryCondtionNumber).xRange = {0,0};
        boundaryConditions->at(boundaryCondtionNumber).yRange = {0,0};
        boundaryConditions->at(boundaryCondtionNumber).zRange = {0,0};

        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "type",       PETSC_INT,    NULL, &boundaryConditions->at(boundaryCondtionNumber).type));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "xRangeMin",  PETSC_SCALAR, NULL, &(boundaryConditions->at(boundaryCondtionNumber).xRange[0]) ));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "xRangeMax",  PETSC_SCALAR, NULL, &(boundaryConditions->at(boundaryCondtionNumber).xRange[1]) ));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "yRangeMin",  PETSC_SCALAR, NULL, &(boundaryConditions->at(boundaryCondtionNumber).yRange[0]) ));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "yRangeMax",  PETSC_SCALAR, NULL, &(boundaryConditions->at(boundaryCondtionNumber).yRange[1]) ));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "zRangeMin",  PETSC_SCALAR, NULL, &(boundaryConditions->at(boundaryCondtionNumber).zRange[0]) ));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "zRangeMax",  PETSC_SCALAR, NULL, &(boundaryConditions->at(boundaryCondtionNumber).zRange[1]) ));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "value",      PETSC_SCALAR, NULL, &boundaryConditions->at(boundaryCondtionNumber).value));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "dofX",       PETSC_BOOL,   NULL, &dofX));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "dofY",       PETSC_BOOL,   NULL, &dofY));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "dofZ",       PETSC_BOOL,   NULL, &dofZ));

        if (dofX)
            boundaryConditions->at(boundaryCondtionNumber).degreesOfFreedom.insert(0);
        if (dofY)
            boundaryConditions->at(boundaryCondtionNumber).degreesOfFreedom.insert(1);
        if (dofZ)
            boundaryConditions->at(boundaryCondtionNumber).degreesOfFreedom.insert(2);

        boundaryCondtionNumber++;
        parent = settingsGroup + "/boundaryConditions/boundaryCondition_" + std::to_string(boundaryCondtionNumber);
        PetscCall(PetscViewerHDF5HasGroup(saveFileHDF5, parent.c_str(), &hasDataset));
    }

    return 0;
}

PetscErrorCode FileManager::saveDomainInformation(PetscViewer saveFileHDF5, DomainCoordinates domain)
{
    std::string parent = settingsGroup + "/domain";
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "xMinimum", PETSC_SCALAR, &domain.xMinimum));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "xMaximum", PETSC_SCALAR, &domain.xMaximum));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "yMinimum", PETSC_SCALAR, &domain.yMinimum));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "yMaximum", PETSC_SCALAR, &domain.yMaximum));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "zMinimum", PETSC_SCALAR, &domain.zMinimum));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, parent.c_str(), "zMaximum", PETSC_SCALAR, &domain.zMaximum));

    return 0;
}

PetscErrorCode FileManager::getDomainInformation(PetscViewer saveFileHDF5, DomainCoordinates *domain)
{
    std::string parent = settingsGroup + "/domain";
    PetscBool hasDataset = PETSC_FALSE;
    PetscCall(PetscViewerHDF5HasGroup(saveFileHDF5, parent.c_str(), &hasDataset));

    if (hasDataset)
    {
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "xMinimum", PETSC_SCALAR, NULL, &(domain->xMinimum)));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "xMaximum", PETSC_SCALAR, NULL, &(domain->xMaximum)));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "yMinimum", PETSC_SCALAR, NULL, &(domain->yMinimum)));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "yMaximum", PETSC_SCALAR, NULL, &(domain->yMaximum)));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "zMinimum", PETSC_SCALAR, NULL, &(domain->zMinimum)));
        PetscCall(PetscViewerHDF5ReadAttribute(saveFileHDF5, parent.c_str(), "zMaximum", PETSC_SCALAR, NULL, &(domain->zMaximum)));
    }

    return 0;
}
