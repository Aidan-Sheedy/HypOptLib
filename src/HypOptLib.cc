

#include "HypOptLib.h"
#include "Hyperoptimization.h"
#include "SensitivitiesWrapper.h"
#include "FilterWrapper.h"
#include "LagrangeMultiplier.h"
#include "FileManager.h"

#include <petscsystypes.h>
#include <random>

void HypOptLib::init(   bool                    randomStartingValues,
                        bool                    saveHamiltonian,
                        double                  initialTemperature,
                        double                  penalty,
                        double                  minimumFilterRadius,
                        double                  volumeFraction,
                        double                  timestep,
                        uint32_t                noseHooverChainOrder,
                        uint32_t                maximumIterations,
                        std::vector<double>    *iterationSaveRange,
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

    this->opt = new TopOpt(xMeshDimension, yMeshDimension, zMeshDimension, penalty, minimumFilterRadius);

    // STEP 2: THE PHYSICS
    this->physics = new LinearElasticity(opt->da_nodes);

    // STEP 3: THE FILTERING
    this->filter = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);

    this->saveHamiltonian       = saveHamiltonian;
    this->initialTemperature    = initialTemperature;
    this->timestep              = timestep;
    this->noseHooverChainOrder  = noseHooverChainOrder;
    this->maximumIterations     = maximumIterations;
    this->iterationSaveRange    = *iterationSaveRange;
    this->volumeFraction        = volumeFraction;
    this->gridDimensions        = *gridDimensions;

    PetscErrorCode status = initNoRestart(randomStartingValues);
}

PetscErrorCode HypOptLib::initNoRestart(bool randomStartingValues)
{
    PetscPrintf(PETSC_COMM_WORLD, "\n\n######################## Initialization ########################\n");

    /* Clone Petsc settings from topopt (it initializes everything) */
    PetscCall(VecDuplicate(opt->x, &initialPositions));
    PetscCall(VecDuplicate(opt->x, &initialVelocities));

    if (randomStartingValues)
    {
        /* Get underlying vector array of initial positions */
        PetscScalar *initialValues;
        PetscInt     localSize;
        PetscCall(VecGetArray(initialPositions, &initialValues));
        PetscCall(VecGetLocalSize(initialPositions, &localSize));

        PetscPrintf(PETSC_COMM_WORLD, "# Initializing random number genrator.\n");

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
        
        PetscPrintf(PETSC_COMM_WORLD, "# Assigning Random Numbers.\n");
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

    PetscCall(VecSet(initialVelocities, std::sqrt(initialTemperature)));

    PetscPrintf(PETSC_COMM_WORLD, "# Initialized.\n\n");

    return 0;
}

void HypOptLib::startLoop()
{
    SensitivitiesWrapper    sensitivities(physics, opt);
    FilterWrapper           wrappedFilter(filter, opt->m);
    LagrangeMultiplier      lagrangianMultiplier(filter, opt);
    FileManager output;

    PetscPrintf(PETSC_COMM_WORLD, "######################## Readying Design Loop ########################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Initialing File Manager\n");

    output.initializeHDF5(volumeFraction,
                          timestep,
                          initialTemperature,
                          gridDimensions[0],
                          gridDimensions[1],
                          gridDimensions[2],
                          opt->penal,
                          opt->rmin,
                          maximumIterations,
                          maximumIterations,
                          noseHooverChainOrder);

    PetscPrintf(PETSC_COMM_WORLD, "# Initialing Solver\n");

    Hyperoptimization solver;
    PetscErrorCode status = solver.init(&sensitivities,
                                        &wrappedFilter,
                                        lagrangianMultiplier,
                                        initialTemperature,
                                        initialPositions,
                                        initialVelocities,
                                        noseHooverChainOrder,
                                        maximumIterations,
                                        timestep,
                                        &output,
                                        iterationSaveRange[1],
                                        saveHamiltonian);

    if (0 != status)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Failed to initialize solver, aborting with error code %i\n", status);
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "# Initialized, starting design loop\n");

        double t1 = MPI_Wtime();
        status = solver.runDesignLoop();
        double t2 = MPI_Wtime();

        if (0 != status)
        {
            PetscPrintf(PETSC_COMM_WORLD, "Failed to run solver, aborting with error code %i\n", status);
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD, "# Done loop!\nTotal Runtime: %f\nAverage Iteration Time: %f\n***Note these timings include file saving***", t2 - t1, (t2-t1)/opt->maxItr);
            PetscPrintf(PETSC_COMM_WORLD, "\n\n# Hyperoptimization complete. Cleaning up resources. This HypOptLib object is no longer accessible.\n");
            PetscPrintf(PETSC_COMM_WORLD, "###################################################################\n");
        }
    }

    delete filter;
    delete opt;
    delete physics;
    PetscFinalize();
}