/**
 * @file Hyperoptimization.h
 *
 * This file describes the layout of the Hyperoptimization class.
 *
 * @todo LICENSE
**/

#pragma once

#include <petsc.h>
#include <vector>
#include <string>

#include "LinearElasticity.h"
#include "HypOptParameters.h"

class FileManager
{
    public:
        // FileManager(){}

        FileManager(LinearElasticity *physics) : physics(physics) {}

        // initializeHDF5(std::string fileName, std::vector<std::string> settings);

        static PetscErrorCode HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<PetscScalar> vector, const char * vectorName);

        static PetscErrorCode HDF5GetSavedVec(std::string filePath, std::string location, Vec vector);

        static PetscErrorCode getHDF5Settings(  std::string             filePath,
                                                PetscInt               *numSavedIterations,
                                                PetscInt               *noseHooverChainOrder,
                                                PetscScalar            *volumeFraction,
                                                PetscScalar            *timestep,
                                                PetscScalar            *targetTemperature,
                                                PetscScalar            *penalty,
                                                PetscScalar            *minimumFilterRadius,
                                                std::vector<uint32_t>  *gridDimensions);

        PetscErrorCode initializeHDF5(PetscScalar   volfrac,
                                      PetscScalar   timestep,
                                      PetscScalar   temperature,
                                      PetscInt      numberElementsX,
                                      PetscInt      numberElementsY,
                                      PetscInt      numberElementsZ,
                                      PetscScalar   penalty,
                                      PetscScalar   filterRadius,
                                      PetscInt      numberSteps,
                                      PetscInt      numberSamples,
                                      PetscInt      NoseHooverChainOrder,
                                      std::string   filePath);

        static PetscErrorCode getHDF5Vectors(std::string        filePath,
                                            HypOptParameters    finalState,
                                            Vec                 finalStateField);

        static bool doesFileExist(std::string filepath);

        PetscErrorCode saveIteration(PetscInt iteration, Vec positions);

        PetscErrorCode saveFinalState(  bool saveHamiltionan,
                                        HypOptParameters finalState,
                                        std::vector<PetscScalar> hamiltonians,
                                        std::vector<PetscScalar> compliance,
                                        std::vector<PetscScalar> temperatures,
                                        std::vector<PetscScalar> LagrangeMultipliers,
                                        std::vector<PetscScalar> iterationTimes);

        std::string getSaveFilePath();
        std::string getDataGroup();

    private:
        std::string autoAppendFilePath(std::string fileName, std::string fileExtension);

        // const std::string stateGroup = "/Dataset/State";
        std::string saveFilePath;

        const std::string stateGroup = "/Dataset/State";

        const std::string dataGroup = "/Dataset";

        LinearElasticity *physics;
};
