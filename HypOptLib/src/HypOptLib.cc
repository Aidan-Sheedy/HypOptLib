/***************************************************************************//**
 * @file HypOptLib.cc
 *
 * Contains the main HypOptLib class, implementation
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
 ******************************************************************************/

#include "HypOptLib.h"
#include "SensitivitiesWrapper.h"
#include "FilterWrapper.h"
#include "LagrangeMultiplier.h"
#include "FileManager.h"
#include "HypOptParameters.h"

#include <petscsystypes.h>
#include <random>
#include <fstream>
#include <iostream>

uint32_t HypOptLib::newRun( std::vector<uint32_t>  *iterationSaveRange,
                            std::vector<uint32_t>  *gridDimensions)
{
    /* Check for valid input parameters */
    if (3 != gridDimensions->size())
    {
        throw HypOptException("Invalid grid dimensions. Must provide a list of exactly three dimensions; x, y, and z.");
    }

    if (2 != iterationSaveRange->size())
    {
        throw HypOptException("Invalid iteration save range. Must be a list of two values, first smaller than the last.");
    }

    // Initialize PETSc / MPI and pass input arguments to PETSc
    PetscInitializeNoArguments();

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)

    PetscInt xMeshDimension = gridDimensions->at(0)+1;
    PetscInt yMeshDimension = gridDimensions->at(1)+1;
    PetscInt zMeshDimension = gridDimensions->at(2)+1;

    TopOpt* opt = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, penalty, minimumFilterRadius);

    // STEP 2: THE PHYSICS
    LinearElasticity* physics = new LinearElasticity(opt->da_nodes);

    // STEP 3: THE FILTERING
    Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);

    PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "########################### Hyperoptimization ##########################\n");

    /* Clone Petsc settings from topopt (it initializes everything) */
    Vec initialPositions;
    Vec initialVelocities;
    PetscCall(VecDuplicate(opt->x, &initialPositions));
    PetscCall(VecDuplicate(opt->x, &initialVelocities));

    if (randomStartingValues)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Randomizing initial values.\n");
        randomizeStartingVectors(initialPositions, initialVelocities);
    }
    else if (initialConditionsFromFile)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Loading initial values from file: %s\n", initialConditionsFile.c_str());
        FileManager::loadInitialConditions(initialPositions, initialVelocities, initialConditionsFile);
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Using all identical initial conditions.\n");
        PetscCall(VecSet(initialPositions, volumeFraction));
        PetscCall(VecSet(initialVelocities, std::sqrt(targetTemperature)));
    }

    SensitivitiesWrapper    sensitivities(physics, opt->Emin, opt->Emax, penalty, volumeFraction);
    FilterWrapper           wrappedFilter(filter, opt->m);
    LagrangeMultiplier      lagrangianMultiplier(wrappedFilter, volumeFraction);
    FileManager             output(physics);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing File Manager\n");

    output.initializeHDF5(volumeFraction,
                          timestep,
                          targetTemperature,
                          gridDimensions->at(0),
                          gridDimensions->at(1),
                          gridDimensions->at(2),
                          opt->penal,
                          opt->rmin,
                          maximumIterations,
                          maximumIterations,
                          noseHooverChainOrder,
                          this->savePath,
                          randomStartingValues,
                          initialConditionsFile);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    Hyperoptimization solver;
    PetscErrorCode status = solver.init(sensitivities,
                                        wrappedFilter,
                                        lagrangianMultiplier,
                                        targetTemperature,
                                        initialPositions,
                                        initialVelocities,
                                        noseHooverChainOrder,
                                        maximumIterations,
                                        timestep,
                                        &output,
                                       *iterationSaveRange,
                                        saveHamiltonian,
                                        volumeFraction,
                                        saveFrequency,
                                        maxSimTime);

    if (0 != status)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Failed to initialize solver, aborting with error code %i\n", status);
    }
    else
    {
        runLoop(solver, maximumIterations, output);
    }

    PetscCall(VecDestroy(&initialPositions));
    PetscCall(VecDestroy(&initialVelocities));

    delete filter;
    delete opt;
    delete physics;
    PetscFinalize();

    return 0;
}

uint32_t HypOptLib::restartRun( std::string restartPath,
                                std::vector<uint32_t> *iterationSaveRange)
{
    /* Try and find the file to see if it exists */
    if (!FileManager::doesFileExist(restartPath))
    {
        throw HypOptException("Failed to find desired restart file.");
    }

    if (2 != iterationSaveRange->size())
    {
        throw HypOptException("Invalid iteration save range. Must be a list of two values, first smaller than the last.");
    }

    // Initialize PETSc / MPI and pass input arguments to PETSc
    PetscInitializeNoArguments();

    PetscInt    noseHooverChainOrder;
    PetscScalar volumeFraction;
    PetscScalar timestep;
    PetscScalar targetTemperature;
    PetscScalar penalty;
    PetscScalar minimumFilterRadius;
    std::vector<uint32_t>   gridDimensions;

    FileManager::getHDF5Settings(restartPath,
                                &noseHooverChainOrder,
                                &volumeFraction,
                                &timestep,
                                &targetTemperature,
                                &penalty,
                                &minimumFilterRadius,
                                &gridDimensions);

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)

    PetscInt xMeshDimension = gridDimensions.at(0)+1;
    PetscInt yMeshDimension = gridDimensions.at(1)+1;
    PetscInt zMeshDimension = gridDimensions.at(2)+1;

    TopOpt* opt                 = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, penalty, minimumFilterRadius);
    LinearElasticity* physics   = new LinearElasticity(opt->da_nodes);
    Filter* filter              = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);

    PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "########################### Hyperoptimization ##########################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Restarting with file: %s\n", restartPath.c_str());

    HypOptParameters finalState;
    Vec finalStateField;

    physics->DuplicateStateField(&finalStateField);

    PetscCall(VecDuplicate(opt->x, &(finalState.position)));
    PetscCall(VecDuplicate(opt->x, &(finalState.velocity)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.evenNoseHooverPosition)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.evenNoseHooverVelocity)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.oddNoseHooverVelocity)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.oddNoseHooverPosition)));

    FileManager::getFinalStateVectors(restartPath,
                                finalState,
                                finalStateField);

    physics->CopyVecToStateField(finalStateField);

    SensitivitiesWrapper    sensitivities(physics, opt->Emin, opt->Emax, penalty, volumeFraction);
    FilterWrapper           wrappedFilter(filter, opt->m);
    LagrangeMultiplier      lagrangianMultiplier(wrappedFilter, volumeFraction);
    FileManager             output(physics);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing File Manager\n");

    output.initializeHDF5(volumeFraction,
                          timestep,
                          targetTemperature,
                          gridDimensions[0],
                          gridDimensions[1],
                          gridDimensions[2],
                          opt->penal,
                          opt->rmin,
                          maximumIterations,
                          maximumIterations,
                          noseHooverChainOrder,
                          this->savePath,
                          randomStartingValues,
                          initialConditionsFile);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    Hyperoptimization solver;
    PetscErrorCode status = solver.init(sensitivities,
                                        lagrangianMultiplier,
                                        wrappedFilter,
                                       &output,
                                        noseHooverChainOrder,
                                        timestep,
                                        targetTemperature,
                                        maximumIterations,
                                       *iterationSaveRange,
                                        finalState.position,
                                        finalState.velocity,
                                        finalState.evenNoseHooverPosition,
                                        finalState.evenNoseHooverVelocity,
                                        finalState.oddNoseHooverPosition,
                                        finalState.oddNoseHooverVelocity,
                                        saveHamiltonian,
                                        volumeFraction,
                                        saveFrequency,
                                        maxSimTime);

    if (0 != status)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Failed to initialize solver, aborting with error code %i\n", status);
    }
    else
    {
        runLoop(solver, maximumIterations, output);
    }


    PetscCall(VecDestroy(&(finalState.evenNoseHooverPosition)));
    PetscCall(VecDestroy(&(finalState.evenNoseHooverVelocity)));
    PetscCall(VecDestroy(&(finalState.oddNoseHooverPosition)));
    PetscCall(VecDestroy(&(finalState.oddNoseHooverVelocity)));
    PetscCall(VecDestroy(&(finalState.position)));
    PetscCall(VecDestroy(&(finalState.velocity)));
    PetscCall(VecDestroy(&finalStateField));

    delete filter;
    delete opt;
    delete physics;
    PetscFinalize();

   return 0;
}

void HypOptLib::generateRandomInitialConditionsFile(std::vector<uint32_t> *gridDimensions, std::string filePath)
{
    /* Check for valid input parameters */
    if (3 != gridDimensions->size())
    {
        throw HypOptException("Invalid grid dimensions. Must provide a list of exactly three dimensions; x, y, and z.");
    }

    PetscErrorCode status = 0;

    /* Initialize PETSc / MPI and pass input arguments to PETSc */
    PetscInitializeNoArguments();

    /* Initialize TopOpt object */
    PetscInt xMeshDimension = gridDimensions->at(0)+1;
    PetscInt yMeshDimension = gridDimensions->at(1)+1;
    PetscInt zMeshDimension = gridDimensions->at(2)+1;
    TopOpt* opt = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, penalty, minimumFilterRadius);

    /* Clone Petsc vectors from topopt */
    Vec initialPositions;
    Vec initialVelocities;
    VecDuplicate(opt->x, &initialPositions);
    VecDuplicate(opt->x, &initialVelocities);

    PetscPrintf(PETSC_COMM_WORLD, "########################### Hyperoptimization ##########################\n");

    PetscPrintf(PETSC_COMM_WORLD, "# Generating Random Starting Vectors given:\n");
    PetscPrintf(PETSC_COMM_WORLD, "# - Temperature: %f\n", targetTemperature);
    PetscPrintf(PETSC_COMM_WORLD, "# - Dimensions: (%d, %d, %d)\n", gridDimensions->at(0), gridDimensions->at(1), gridDimensions->at(2));

    status = randomizeStartingVectors(initialPositions, initialVelocities);

    if (0 != status)
    {
        throw HypOptException((std::string("Error randomizing starting values: ") + std::to_string(status)).c_str());
    }

    PetscPrintf(PETSC_COMM_WORLD, "# Generated!\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Saving to file: %s\n", filePath.c_str());

    FileManager::saveInitialConditions(initialPositions, initialVelocities, filePath);

    PetscPrintf(PETSC_COMM_WORLD, "# Done!\n");
    PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");

    VecDestroy(&initialPositions);
    VecDestroy(&initialVelocities);
    delete opt;
}

PetscErrorCode HypOptLib::randomizeStartingVectors(Vec initialPosition, Vec initialVelocity)
{
    /* Get underlying vector array of initial positions */
    PetscScalar *initialValues;
    PetscInt     localSize;
    PetscInt     globalSize;
    PetscCall(VecGetArray(initialPosition, &initialValues));
    PetscCall(VecGetLocalSize(initialPosition, &localSize));
    PetscCall(VecGetSize(initialPosition, &globalSize));

    /* Initialize random number generator */
    std::default_random_engine generator;
    std::uniform_real_distribution<PetscScalar> distribution;
    if (0.5 < volumeFraction)
    {
        distribution = std::uniform_real_distribution<PetscScalar>(2 * volumeFraction - 1, 1);
    }
    else
    {
        distribution = std::uniform_real_distribution<PetscScalar>(0, 2 * volumeFraction);
    }

    /* Assign uniform random distribution as initial values */
    for (PetscInt i = 0; i < localSize; i++)
    {
        initialValues[i] = distribution(generator); //volumeFraction;
    }

    /* Restore Vector */
    PetscCall(VecRestoreArray(initialPosition, &initialValues));
    PetscCall(VecSet(initialVelocity, std::sqrt(targetTemperature)));

    /* Assign velocity initial as Maxwell Boltzmann distribution */
    PetscCall(VecGetArray(initialVelocity, &initialValues));
    PetscCall(VecGetLocalSize(initialVelocity, &localSize));

    std::default_random_engine generator2;
    std::uniform_real_distribution<PetscScalar> distribution2;

    // distribution2.reset();
    distribution2 = std::uniform_real_distribution<PetscScalar>(-0.5, 0.5);
    for (PetscInt i = 0; i < localSize; i++)
    {
        initialValues[i] = distribution2(generator); //volumeFraction;
    }

    PetscCall(VecRestoreArray(initialVelocity, &initialValues));

    Vec temp;
    PetscScalar sumV;
    PetscScalar sumV2;

    PetscCall(VecDuplicate(initialVelocity, &temp));
    PetscCall(VecCopy(initialVelocity, temp));
    PetscCall(VecPointwiseMult(temp, initialVelocity, initialVelocity));

    PetscCall(VecMean(initialVelocity, &sumV));
    PetscCall(VecMean(temp, &sumV2));

    PetscScalar offset = sqrt(targetTemperature / sumV2);

    PetscCall(VecShift(initialVelocity, -sumV));
    PetscCall(VecScale(initialVelocity, offset));

    VecDestroy(&temp);

    return 0;
}

PetscErrorCode HypOptLib::runLoop(Hyperoptimization solver, PetscInt numItr, FileManager output)
{
    PetscPrintf(PETSC_COMM_WORLD, "# Initialized, starting design loop\n");

    if (variableTimestep)
    {
        PetscPrintf(PETSC_COMM_WORLD, "#   - Enabling variable timestep\n");
        solver.enableVariableTimestep(timestepConstantAlpha, timestepConstantBeta, diffusionConstant);
    }

    double t1 = MPI_Wtime();
    PetscErrorCode status = solver.runDesignLoop();
    double t2 = MPI_Wtime();

    if (0 != status)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Failed to run solver, aborting with error code %i\n", status);
    }
    else
    {
        output.saveFinalState(  solver.getSaveHamiltonian(),
                                solver.getFinalState(),
                                solver.getHamiltonians(),
                                solver.getCompliance(),
                                solver.getTemperatures(),
                                solver.getLagrangeMultipliers(),
                                solver.getIterationTimes(),
                                solver.getTimesteps(),
                                solver.getEnergyErrors(),
                                solver.getVolFracs());

        PetscPrintf(PETSC_COMM_WORLD, "# Total Runtime: %f\n# Average Iteration Time: %f\n# ***Note these timings include file saving***", t2 - t1, (t2-t1)/numItr);
        PetscPrintf(PETSC_COMM_WORLD, "\n# Hyperoptimization complete, cleaning up resources.\n");
        PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");
    }

   return 0;
}