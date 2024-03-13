/***************************************************************************//**
 * @file HypOptLib.cc
 *
 * Contains the main HypOptLib class, implementation
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

#include "HypOptLib.h"
#include "SensitivitiesWrapper.h"
#include "FilterWrapper.h"
#include "LagrangeMultiplier.h"
#include "FileManager.h"

#include <petscsystypes.h>
#include <random>
#include <fstream>
#include <iostream>

uint32_t HypOptLib::newRun(std::vector<uint32_t>  *iterationSaveRange)
{
    /* Check for valid input parameters */
    if (2 != iterationSaveRange->size())
    {
        throw HypOptException("Invalid iteration save range. Must be a list of two values, first smaller than the last.");
    }

    if ( (0 == domain.xMinimum) &&
         (0 == domain.xMaximum) &&
         (0 == domain.yMinimum) &&
         (0 == domain.yMaximum) &&
         (0 == domain.zMinimum) &&
         (0 == domain.zMaximum) )
    {
        throw HypOptException("Grid topology not set. Call HypOptLib::setGridProperties with valid properties.");
    }

    if (0 == boundaryConditions.size() || 3 != gridDimensions.size())
    {
        throw HypOptException("Boundary conditions not set. Please call HypOptLib::setBoundaryConditions with valid boundary conditions.");
    }

    // Initialize PETSc / MPI and pass input arguments to PETSc
    PetscInitializeNoArguments();

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)

    PetscInt xMeshDimension = gridDimensions.at(0)+1;
    PetscInt yMeshDimension = gridDimensions.at(1)+1;
    PetscInt zMeshDimension = gridDimensions.at(2)+1;

    TopOpt* opt = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, multigridLevels, domain);

    // STEP 2: THE PHYSICS
    LinearElasticity* physics = new LinearElasticity(opt->da_nodes, boundaryConditions, maxFeaIterations, multigridLevels);

    // STEP 3: THE FILTERING
    /* Always use filter type 1 -- none others are supported by Hyperoptimization. */
    Filter* filter = new Filter(opt->da_nodes, opt->xPhys, 1, minimumFilterRadius);

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

    SensitivitiesWrapper    sensitivities(physics, this->Emin, this->Emax, penalty, volumeFraction);
    FilterWrapper           wrappedFilter(filter, numConstraints);
    LagrangeMultiplier      lagrangianMultiplier(wrappedFilter, volumeFraction);
    FileManager             output(physics);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing File Manager\n");

    output.initializeHDF5(volumeFraction,
                          timestep,
                          targetTemperature,
                          gridDimensions.at(0),
                          gridDimensions.at(1),
                          gridDimensions.at(2),
                          penalty,
                          minimumFilterRadius,
                          maximumIterations,
                          saveFrequency,
                          noseHooverChainOrder,
                          this->savePath,
                          randomStartingValues,
                          initialConditionsFile,
                          domain,
                          boundaryConditions);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    printOptions(*iterationSaveRange);

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
                                        printInfo,
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
    DomainCoordinates previousDomain;
    std::vector<BoundaryCondition> prevBoundaryConditions;

    FileManager::getHDF5Settings(restartPath,
                                &noseHooverChainOrder,
                                &volumeFraction,
                                &timestep,
                                &targetTemperature,
                                &penalty,
                                &minimumFilterRadius,
                                &gridDimensions,
                                &previousDomain,
                                &prevBoundaryConditions);

    if (overrideRestartTimestep)
    {
        timestep = this->timestep;
    }

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)

    PetscInt xMeshDimension = gridDimensions.at(0)+1;
    PetscInt yMeshDimension = gridDimensions.at(1)+1;
    PetscInt zMeshDimension = gridDimensions.at(2)+1;

    DomainCoordinates               newDomain;
    std::vector<BoundaryCondition>  newBCs;

    if (restartUseNewMeshes)
    {
        newDomain = this->domain;
        newBCs = this->boundaryConditions;
    }
    else
    {
        newDomain = previousDomain;
        newBCs = prevBoundaryConditions;
    }

    TopOpt* opt                 = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, multigridLevels, newDomain);
    LinearElasticity* physics   = new LinearElasticity(opt->da_nodes, newBCs, maxFeaIterations, multigridLevels);

    /* Always use filter type 1 -- none others are supported by Hyperoptimization. */
    Filter* filter              = new Filter(opt->da_nodes, opt->xPhys, 1, minimumFilterRadius);

    PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "########################### Hyperoptimization ##########################\n");

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

    SensitivitiesWrapper    sensitivities(physics, this->Emin, this->Emax, penalty, volumeFraction);
    FilterWrapper           wrappedFilter(filter, numConstraints);
    LagrangeMultiplier      lagrangianMultiplier(wrappedFilter, volumeFraction);
    FileManager             output(physics);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing File Manager\n");

    output.initializeHDF5(volumeFraction,
                          timestep,
                          targetTemperature,
                          gridDimensions[0],
                          gridDimensions[1],
                          gridDimensions[2],
                          penalty,
                          minimumFilterRadius,
                          maximumIterations,
                          saveFrequency,
                          noseHooverChainOrder,
                          this->savePath,
                          randomStartingValues,
                          initialConditionsFile,
                          newDomain,
                          newBCs);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    printOptions(*iterationSaveRange, "restart", restartPath);

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
                                        printInfo,
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
    DomainCoordinates domain;
    domain.xMinimum = 0;
    domain.xMaximum = 0;
    domain.yMinimum = 0;
    domain.yMaximum = 0;
    domain.zMinimum = 0;
    domain.zMaximum = 0;
    TopOpt* opt = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, multigridLevels, domain);

    /* Clone Petsc vectors from topopt */
    Vec initialPositions;
    Vec initialVelocities;
    VecDuplicate(opt->x, &initialPositions);
    VecDuplicate(opt->x, &initialVelocities);

    PetscPrintf(PETSC_COMM_WORLD, "###################### Generate Initial Conditions #####################\n");

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
        PetscPrintf(PETSC_COMM_WORLD, "# Enabling variable timestep\n");
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
        PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");

        output.saveFinalState(  solver.getSaveHamiltonian(),
                                solver.getFinalState(),
                                solver.getHamiltonians(),
                                solver.getCompliance(),
                                solver.getTemperatures(),
                                solver.getLagrangeMultipliers(),
                                solver.getIterationTimes(),
                                solver.getTimesteps(),
                                solver.getEnergyErrors(),
                                solver.getVolFracs(),
                                solver.getSolverIterationsSensitivity(),
                                solver.getSolverIterationsHamiltonian());

        PetscPrintf(PETSC_COMM_WORLD, "# Total Runtime: %f\n# Average Iteration Time: %f\n# ***Note these timings include file saving***", t2 - t1, (t2-t1)/numItr);
        PetscPrintf(PETSC_COMM_WORLD, "\n# Hyperoptimization complete, cleaning up resources.\n");
        PetscPrintf(PETSC_COMM_WORLD, "########################################################################\n");
    }

   return 0;
}

void HypOptLib::setGridProperties(std::vector<uint32_t> *gridDimensions, DomainCoordinates domain)
{
    /* Check for valid input parameters */
    if (3 != gridDimensions->size())
    {
        throw HypOptException("Invalid grid dimensions. Must provide a list of exactly three dimensions; x, y, and z.");
    }

    if (domain.xMinimum > domain.xMaximum)
    {
        throw HypOptException("Invalid domain x coordinates.");
    }

    if (domain.yMinimum > domain.yMaximum)
    {
        throw HypOptException("Invalid domain y coordinates.");
    }

    if (domain.zMinimum > domain.zMaximum)
    {
        throw HypOptException("Invalid domain z coordinates.");
    }

    /* Make deep copy to avoid Python cleaning up the variable */
    for (auto gridDimension : *gridDimensions)
        this->gridDimensions.push_back(gridDimension);

    this->domain.xMinimum = domain.xMinimum;
    this->domain.xMaximum = domain.xMaximum;
    this->domain.yMinimum = domain.yMinimum;
    this->domain.yMaximum = domain.yMaximum;
    this->domain.zMinimum = domain.zMinimum;
    this->domain.zMaximum = domain.zMaximum;

    return;
}

void HypOptLib::setBoundaryConditions(std::vector<BoundaryCondition> *boundaryConditions)
{
    /* Sanitize boundary conditions */
    if (0 == boundaryConditions->size())
    {
        throw HypOptException("ERROR! No boundary conditions supplied.");
    }

    for (auto boundaryCondition : *boundaryConditions)
    {
        if (boundaryCondition.xRange.size() > 0 && boundaryCondition.xRange.size() != 2)
        {
            throw HypOptException("ERROR! Bad boundary conditions supplied.");
        }
        else if (boundaryCondition.xRange[0] > boundaryCondition.xRange[1])
        {
            throw HypOptException("ERROR! Bad boundary conditions supplied, range minimum must be smaller or equal to the max.");
        }

        if (boundaryCondition.yRange.size() > 0 && boundaryCondition.yRange.size() != 2)
        {
            throw HypOptException("ERROR! Bad boundary conditions supplied.");
        }
        else if (boundaryCondition.yRange[0] > boundaryCondition.yRange[1])
        {
            throw HypOptException("ERROR! Bad boundary conditions supplied, range minimum must be smaller or equal to the max.");
        }

        if (boundaryCondition.zRange.size() > 0 && boundaryCondition.zRange.size() != 2)
        {
            throw HypOptException("ERROR! Bad boundary conditions supplied.");
        }
        else if (boundaryCondition.zRange[0] > boundaryCondition.zRange[1])
        {
            throw HypOptException("ERROR! Bad boundary conditions supplied, range minimum must be smaller or equal to the max.");
        }

        if ( (0 == boundaryCondition.xRange.size()) && 
             (0 == boundaryCondition.yRange.size()) &&
             (0 == boundaryCondition.zRange.size()) )
        {
            throw HypOptException("ERROR! Empty boundary condition supplied.");
        }

        if ( (boundaryCondition.degreesOfFreedom.end() != boundaryCondition.degreesOfFreedom.upper_bound(3)) ||
             (boundaryCondition.degreesOfFreedom.empty())                                                    )
        {
            throw HypOptException("ERROR! Degrees of freedom is not valid.");
        }

        this->boundaryConditions.push_back(boundaryCondition);
    }
}

void HypOptLib::printOptions(std::vector<uint32_t>  iterationSaveRange, std::string options, std::string restartPath)
{

    if (DEBUG <= printInfo)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# - Boundary Conditions:\n");

        for (auto boundaryCondition : boundaryConditions)
        {
            PetscPrintf(PETSC_COMM_WORLD, "#       type:    %i\n", boundaryCondition.type);
            PetscPrintf(PETSC_COMM_WORLD, "#       xRange:  (%.1f,%.1f)\n", boundaryCondition.xRange[0], boundaryCondition.xRange[1]);
            PetscPrintf(PETSC_COMM_WORLD, "#       yRange:  (%.1f,%.1f)\n", boundaryCondition.yRange[0], boundaryCondition.yRange[1]);
            PetscPrintf(PETSC_COMM_WORLD, "#       zRange:  (%.1f,%.1f)\n", boundaryCondition.zRange[0], boundaryCondition.zRange[1]);
            PetscPrintf(PETSC_COMM_WORLD, "#       DOF:     (");
            for (PetscInt dof : boundaryCondition.degreesOfFreedom)
                PetscPrintf(PETSC_COMM_WORLD, "%i,", dof);
            PetscPrintf(PETSC_COMM_WORLD, ")\n#       value:    %.2e\n", boundaryCondition.value);
        }
    }
    if ("restart" == options)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Restarting with file: %s\n", restartPath.c_str());
    }
    PetscPrintf(PETSC_COMM_WORLD, "# ------------ Hyperoptimization Settings ------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "# - targetTemperature:         %.2e\n", targetTemperature);
    if (INFO <= printInfo)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# - noseHooverChainOrder:      %d\n", noseHooverChainOrder);
        PetscPrintf(PETSC_COMM_WORLD, "# - randomStartingValues:      %s\n", (0 == randomStartingValues)        ? "False" : "True");
        PetscPrintf(PETSC_COMM_WORLD, "# - saveHamiltonian:           %s\n", (0 == saveHamiltonian)             ? "False" : "True");
    }
    PetscPrintf(PETSC_COMM_WORLD, "# --------------- Simulation Settings ----------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "# - timestep:                  %.2e\n", timestep);
    PetscPrintf(PETSC_COMM_WORLD, "# - maximumIterations:         %d\n", maximumIterations);
    if (std::numeric_limits<double>::max() > maxSimTime)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# - maxSimTime:                %.2e\n", maxSimTime);
    }
    if (variableTimestep && (DEBUG <= printInfo))
    {
        PetscPrintf(PETSC_COMM_WORLD, "# - timestepConstantAlpha:     %.2e\n", timestepConstantAlpha);
        PetscPrintf(PETSC_COMM_WORLD, "# - timestepConstantBeta:      %.2e\n", timestepConstantBeta);
        PetscPrintf(PETSC_COMM_WORLD, "# - diffusionConstant:         %.2e\n", diffusionConstant);
    }
    PetscPrintf(PETSC_COMM_WORLD, "# - iterationSaveRange:        (%i, %i)\n", iterationSaveRange[0], iterationSaveRange[1]);
    PetscPrintf(PETSC_COMM_WORLD, "# - saveFrequency:             %d\n", saveFrequency);
    if (DEBUG <= printInfo)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# - restartUseNewMeshes:       %s\n", (0 == restartUseNewMeshes)         ? "False" : "True");
        PetscPrintf(PETSC_COMM_WORLD, "# - overrideRestartTimestep:   %s\n", (0 == overrideRestartTimestep)     ? "False" : "True");
    }
    PetscPrintf(PETSC_COMM_WORLD, "# - variableTimestep:          %s\n", (0 == variableTimestep)            ? "False" : "True");
    PetscPrintf(PETSC_COMM_WORLD, "# - initialConditionsFromFile: %s\n", (0 == initialConditionsFromFile)   ? "False" : "True");
    PetscPrintf(PETSC_COMM_WORLD, "# - printInfo:                 %d\n", printInfo);

    if (INFO <= printInfo)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# ---------------- Physical Settings -----------------\n");
        PetscPrintf(PETSC_COMM_WORLD, "# - volumeFraction:            %.2e\n", volumeFraction);
        PetscPrintf(PETSC_COMM_WORLD, "# - minimumFilterRadius:       %.2e\n", minimumFilterRadius);
        PetscPrintf(PETSC_COMM_WORLD, "# - penalty:                   %.2e\n", penalty);
        if (DEBUG <= printInfo)
        {
            PetscPrintf(PETSC_COMM_WORLD, "# - Emin:                      %.2e\n", Emin);
            PetscPrintf(PETSC_COMM_WORLD, "# - Emax:                      %.2e\n", Emax);
            PetscPrintf(PETSC_COMM_WORLD, "# - maxFeaIterations:          %d\n", maxFeaIterations);
            PetscPrintf(PETSC_COMM_WORLD, "# - multigridLevels:           %d\n", multigridLevels);
        }
        PetscPrintf(PETSC_COMM_WORLD, "# ----------------------------------------------------\n");
    }
}
