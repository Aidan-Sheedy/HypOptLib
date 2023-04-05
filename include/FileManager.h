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
        FileManager() = delete;

        // initializeHDF5(std::string fileName, std::vector<std::string> settings);

        static PetscErrorCode HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<PetscScalar> vector, const char * vectorName);

        static PetscErrorCode HDF5GetSavedVec(std::string filePath, Vec *vector);


    private:
        // const std::string stateGroup = "/Dataset/State";



};