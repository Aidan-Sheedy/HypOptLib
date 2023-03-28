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



    // private:
        // std::string saveFilePathHDF5;

};