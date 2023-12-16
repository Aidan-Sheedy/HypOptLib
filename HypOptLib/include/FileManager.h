/***************************************************************************//**
 * @file FileManager.h
 *
 * This file describes the layout of the FileManager class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#pragma once

#include <petsc.h>
#include <vector>
#include <string>
#include "LinearElasticity.h"
#include "HypOptParameters.h"

class FileManager
{
    public:

        /** Constructor. Saves pointer to provided LinearElasticity object.
         * 
         * @param physics pointer to the LinearElasticity object that will be used by the hyperoptimization algorithm.
         */
        FileManager(LinearElasticity *physics) : physics(physics) {}

        /**
         * Creates an HDF5 file for iteration with the provided parameters. This must be called before
         * starting the iteration loop.
         * 
         * @param volfrac volume fraction setting
         * @param timestep timestep setting
         * @param temperature system temperature setting
         * @param numberElementsX cells in the x direction
         * @param numberElementsY cells in the y direction
         * @param numberElementsZ cells in the z direction
         * @param penalty optimization penalty power setting
         * @param filterRadius minimum filter radius setting
         * @param numberSteps number of iterations to save
         * @param numberSamples number of samples to save @todo figure out step vs sample
         * @param NoseHooverChainOrder number of particles in the Nose Hoover chain
         * @param filePath path to save the HDF5 output file to
         * 
         * @return 0 or PetscError.
         */
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
                                      std::string   filePath,
                                      bool          randomStartingValues,
                                      std::string   initialConditionsFile);

        /**
         * Saves an iteration of the optimized positions.
         * 
         * @param iteration step to save
         * @param positions grid of scalar position values at the given iteration
         * 
         * @return 0 or PetscError.
         */
        PetscErrorCode saveIteration(PetscInt iteration, Vec positions);

        /**
         * Saves the final state of the system for restarting functionality. Additionally saves metadata from each iteration.
         * 
         * @param saveHamiltionan option to save the Hamiltonian and Compliance metadata.
         * @param finalState final state of the system, including all particle and NH particle positions and velocities.
         * @param hamiltonians Hamiltonian at each iteration (if saveHamiltionan option is false this can be NULL).
         * @param compliance compliance at each iteration (if saveHamiltionan option is false this can be NULL).
         * @param temperatures system temperature at each iteration.
         * @param LagrangeMultipliers Lagrangian Multipliers at each iteration.
         * @param iterationTimes time to solve each iteration.
         * 
         * @return 0 or PetscError.
         */
        PetscErrorCode saveFinalState(  bool                     saveHamiltionan,
                                        HypOptParameters         finalState,
                                        std::vector<PetscScalar> hamiltonians,
                                        std::vector<PetscScalar> compliance,
                                        std::vector<PetscScalar> temperatures,
                                        std::vector<PetscScalar> LagrangeMultipliers,
                                        std::vector<PetscScalar> iterationTimes,
                                        std::vector<PetscScalar> timesteps,
                                        std::vector<PetscScalar> energyErrors,
                                        std::vector<PetscScalar> volFracs);

        /**
         * Save file path accessor.
         * 
         * @returns the path to the HDF5 save file.
         */
        std::string getSaveFilePath();

        /** 
         * Data group accessor.
         * 
         * @returns the the name of the data group.
         */
        std::string getDataGroup();

        /**
         * Helper function to determine if a file exists.
         * 
         * Uses the std::ifstream::open function and checks if the fail opened succesfully.
         * This method is imperfect if the file fails to open for another reason, but should
         * work for most cases without getting too complicated.
         *
         * @param filePath the path to check.
         * 
         * @returns true if the file exists, false otherwise.
         */
        static bool doesFileExist(std::string filePath);

        /**
         * Helper function to get the final state from a given restart file.
         * 
         * @param filePath the restart file to parse.
         * @param finalState [out] structure to fill with with the state parameters (particle positions and velocities)
         * @param finalStateField [out] final state field vector to pass to the Linear Elasticity class.
         * 
         * @returns 0 for success, PetscError otherwise.
         */
        static PetscErrorCode getFinalStateVectors( std::string         filePath,
                                                    HypOptParameters    finalState,
                                                    Vec                 finalStateField);

        /**
         * Helper function which converts a std::vector object into a Petsc vector before saving it to the
         * provided HDF5 viewer.
         * 
         * Currently only needed for PetscScalar objects, though this could easily be overloaded or templated
         * for more support if needed.
         * 
         * @param HDF5saveFile Petsc Viewer set up to the desired file and data group.
         * @param vector scalar std::vector object to save.
         * @param vectorName name under which to save the vector in the HDF5 file.
         * 
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode HDF5SaveStdVector(PetscViewer HDF5saveFile, std::vector<PetscScalar> vector, const char * vectorName);

        /**
         * Helper function to retrieve a vector from a given HDF5 file.
         *
         * @todo this function should probably be overloaded to pass a pre-preppared PetscViewer
         * for efficiency (although this is far from the slowest part of the code).
         * 
         * @param filePath the file from which to acquire the desired vector.
         * @param location HDF5 group where the vector is saved
         * @param vector   [out] the vector to populate
         * 
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode HDF5GetSavedVec(std::string filePath, std::string location, Vec vector);

        /**
         * Helper function to retrieve the simulation settings froma  given restart file.
         * 
         * @param filePath the restart file to parse.
         * @param noseHooverChainOrder [out]
         * @param volumeFraction [out]
         * @param timestep [out]
         * @param targetTemperature [out]
         * @param penalty [out]
         * @param minimumFilterRadius [out]
         * @param gridDimensions [out]
         */
        static PetscErrorCode getHDF5Settings(  std::string             filePath,
                                                PetscInt               *noseHooverChainOrder,
                                                PetscScalar            *volumeFraction,
                                                PetscScalar            *timestep,
                                                PetscScalar            *targetTemperature,
                                                PetscScalar            *penalty,
                                                PetscScalar            *minimumFilterRadius,
                                                std::vector<uint32_t>  *gridDimensions);

        static PetscErrorCode saveInitialConditions(Vec positions, Vec velocities, std::string filePath);

        static PetscErrorCode loadInitialConditions(Vec positions, Vec velocities, std::string filePath);

        PetscErrorCode saveIteration(PetscInt iteration, Vec positions, PetscInt timestep);
        PetscErrorCode finishIterating();
        PetscErrorCode preparetoIterate(Vec *positions, PetscInt firstTimestep);

    private:
        /**
         * Private helper function which finds a unique file name in the given path.
         * 
         * If the given file alreaqdy exists, " (1)" is appended to the name. If it still does 
         * exists, " (1)" is incremented until a file name is found which does not exist.
         * 
         * @todo this could probably be simplified to not require sepparating file path and extension.
         * 
         * @param fileName the file name to append (with no file extension).
         * @param fileExtension the file extension to be used.
         */
        std::string autoAppendFilePath(std::string fileName, std::string fileExtension);

        /* User set attributes */
        std::string saveFilePath;
        LinearElasticity *physics;

        /** File-related constant attributes defined in FileManager.cc **/

        /* HDF5 groups */
        static const std::string stateGroup;
        static const std::string dataGroup;
        static const std::string settingsGroup;

        /* System parameters and settings */
        static const std::string volfracName;
        static const std::string timestepName;
        static const std::string tempName;
        static const std::string xDimName;
        static const std::string yDimName;
        static const std::string zDimName;
        static const std::string penaltyName;
        static const std::string minRadiusName;
        static const std::string numStepName;
        static const std::string numSampleName;
        static const std::string chainOrdName;

        /* Simulation Data */
        static const std::string hamiltonianName;
        static const std::string complainceName;
        static const std::string temperatureName;
        static const std::string lagrangianMultipliarName;
        static const std::string iterationComputeTimeName;
        static const std::string timestepsName;
        static const std::string energyErrorsName;
        static const std::string finalPositionName;
        static const std::string finalVelocityName;
        static const std::string finalEvenNHPositionName;
        static const std::string finalEvenNHVelocityName;
        static const std::string finalOddNHPositionName;
        static const std::string finalOddNHVelocityName;
        static const std::string finalStateFieldName;

        /* Test file names */
        static const std::string initialPositionName;
        static const std::string initialVelocityName;
};
