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

class FileManager
{
    public:
        FileManager(){}

        // initializeHDF5(std::string fileName, std::vector<std::string> settings);

        static PetscErrorCode HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<PetscScalar> vector, const char * vectorName);

        static PetscErrorCode HDF5GetSavedVec(std::string filePath, std::string location, Vec *vector);

        PetscErrorCode initializeHDF5(PetscScalar volfrac,
                                      PetscScalar timestep,
                                      PetscScalar temperature,
                                      PetscInt numberElementsX,
                                      PetscInt numberElementsY,
                                      PetscInt numberElementsZ,
                                      PetscScalar penalty,
                                      PetscScalar filterRadius,
                                      PetscInt numberSteps,
                                      PetscInt numberSamples,
                                      PetscScalar NoseHooverChainOrder);

        PetscErrorCode saveIteration(PetscInt iteration, Vec positions);

        std::string getSaveFilePath();
        std::string getDataGroup();

    private:
        // const std::string stateGroup = "/Dataset/State";
        std::string saveFilePath;

        const std::string stateGroup = "/Dataset/State";

        const std::string dataGroup = "/Dataset";
};
