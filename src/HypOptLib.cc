

#include "HypOptLib.h"
#include "SensitivitiesWrapper.h"
#include "FilterWrapper.h"
#include "LagrangeMultiplier.h"
#include "FileManager.h"
#include "HypOptParameters.h"

#include <petscsystypes.h>
#include <random>
#include <fstream>

uint32_t HypOptLib::newRun( bool                    randomStartingValues,
                            bool                    saveHamiltonian,
                            double                  targetTemperature,
                            double                  penalty,
                            double                  minimumFilterRadius,
                            double                  volumeFraction,
                            double                  timestep,
                            uint32_t                noseHooverChainOrder,
                            uint32_t                maximumIterations,
                            std::vector<uint32_t>  *iterationSaveRange,
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

    PetscPrintf(PETSC_COMM_WORLD, "\n\n######################## Hyperoptimization ########################\n");

    /* Clone Petsc settings from topopt (it initializes everything) */
    Vec initialPositions;
    Vec initialVelocities;
    PetscCall(VecDuplicate(opt->x, &initialPositions));
    PetscCall(VecDuplicate(opt->x, &initialVelocities));

    if (randomStartingValues)
    {
        /* Get underlying vector array of initial positions */
        PetscScalar *initialValues;
        PetscInt     localSize;
        PetscCall(VecGetArray(initialPositions, &initialValues));
        PetscCall(VecGetLocalSize(initialPositions, &localSize));

        PetscPrintf(PETSC_COMM_WORLD, "# Initializing random number genrator\n");

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
        PetscCall(VecRestoreArray(initialPositions, &initialValues));
    }
    else
    {
        PetscCall(VecSet(initialPositions, volumeFraction));
    }

    PetscCall(VecSet(initialVelocities, std::sqrt(targetTemperature)));

    SensitivitiesWrapper    sensitivities(physics, opt);
    FilterWrapper           wrappedFilter(filter, opt->m);
    LagrangeMultiplier      lagrangianMultiplier(filter, opt);
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
                          this->savePath);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    Hyperoptimization solver;
    PetscErrorCode status = solver.init(&sensitivities,
                                        &wrappedFilter,
                                        lagrangianMultiplier,
                                        targetTemperature,
                                        initialPositions,
                                        initialVelocities,
                                        noseHooverChainOrder,
                                        maximumIterations,
                                        timestep,
                                        &output,
                                       *iterationSaveRange,
                                        saveHamiltonian);

    if (0 != status)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Failed to initialize solver, aborting with error code %i\n", status);
    }
    else
    {
        runLoop(solver, maximumIterations, output);
    }

    delete filter;
    delete opt;
    delete physics;
    PetscFinalize();

    return 0;
}

uint32_t HypOptLib::restartRun( std::string filePath,
                            uint32_t maximumIterations,
                            std::vector<uint32_t> *iterationSaveRange,
                            bool saveHamiltonian)
{
    /* Try and find the file to see if it exists */
    if (!FileManager::doesFileExist(filePath))
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
    PetscInt    numSavedIterations;
    PetscScalar volumeFraction;
    PetscScalar timestep;
    PetscScalar targetTemperature;
    PetscScalar penalty;
    PetscScalar minimumFilterRadius;
    std::vector<uint32_t>   gridDimensions;

    FileManager::getHDF5Settings(filePath,
                                &noseHooverChainOrder,
                                &numSavedIterations,
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

    PetscPrintf(PETSC_COMM_WORLD, "\n\n######################## Hyperoptimization ########################\n");

    HypOptParameters finalState;
    Vec finalStateField;

    physics->DuplicateStateField(&finalStateField);

    PetscCall(VecDuplicate(opt->x, &(finalState.position)));
    PetscCall(VecDuplicate(opt->x, &(finalState.velocity)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.evenNoseHooverPosition)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.evenNoseHooverVelocity)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.oddNoseHooverVelocity)));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(finalState.oddNoseHooverPosition)));

    FileManager::getHDF5Vectors(filePath,
                                finalState,
                                finalStateField);

    physics->CopyVecToStateField(finalStateField);

    SensitivitiesWrapper    sensitivities(physics, opt);
    FilterWrapper           wrappedFilter(filter, opt->m);
    LagrangeMultiplier      lagrangianMultiplier(filter, opt);
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
                          this->savePath);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    Hyperoptimization solver;
    PetscErrorCode status = solver.init(&sensitivities,
                                        lagrangianMultiplier,
                                        &wrappedFilter,
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
                                        saveHamiltonian);

    if (0 != status)
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Failed to initialize solver, aborting with error code %i\n", status);
    }
    else
    {
        runLoop(solver, numSavedIterations, output);
    }

    delete filter;
    delete opt;
    delete physics;
    PetscFinalize();

   return 0;
}

PetscErrorCode HypOptLib::runLoop(Hyperoptimization solver, PetscInt numItr, FileManager output)
{
    PetscPrintf(PETSC_COMM_WORLD, "# Initialized, starting design loop\n");

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
                                solver.getIterationTimes());

        PetscPrintf(PETSC_COMM_WORLD, "# Done loop!\n#\n# Total Runtime: %f\n# Average Iteration Time: %f\n# ***Note these timings include file saving***", t2 - t1, (t2-t1)/numItr);
        PetscPrintf(PETSC_COMM_WORLD, "\n# Hyperoptimization complete, cleaning up resources.\n# Call init or restart before using this object again.\n");
        PetscPrintf(PETSC_COMM_WORLD, "###################################################################\n");
    }

   return 0;
}