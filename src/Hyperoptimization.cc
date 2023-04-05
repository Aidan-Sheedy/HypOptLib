/**
 * @file Hyperoptimization.h
 * 
 * This file describes the layout of the Hyperoptimization class.
 * 
 * @todo LICENSE
**/

#include "Hyperoptimization.h"

#include "PetscExtensions.h"
#include "FileManager.h"
#include <petscviewerhdf5.h>

#include <cmath>
#include <string>
#include <filesystem>
#include <ctime>
#include <chrono>

/**
 * @todo ALL CALLS TO VecSetValues MUST be followed up by "VecAssemblyBegin()" and "VecAssemblyEnd()"!!!!
 *       GO THROUGH THE WHOLE CODE AND UPDATE THIS.
 * 
 * @todo ALL CALLS TO VecGetArrayRead MUST be followed up by "VecRestoreArrayRead()"!!!!!!!!!
*/
PetscErrorCode Hyperoptimization::init( LinearElasticity* physics,
                                        TopOpt* opt,
                                        Filter* filter,
                                        LagrangeMultiplier lagMult,
                                        PetscScalar temperature,
                                        Vec initialPositions,
                                        PetscScalar NHChainOrder,
                                        PetscInt numIterations,
                                        PetscScalar timestep)
{
    return this->init(  physics,
                        opt,
                        filter,
                        lagMult,
                        temperature,
                        initialPositions,
                        NHChainOrder,
                        numIterations,
                        timestep,
                        2000,
                        true); /* Save 2000 iterations by defaultDefault */
}

PetscErrorCode Hyperoptimization::init( LinearElasticity* physics,
                                        TopOpt* opt,
                                        Filter* filter,
                                        LagrangeMultiplier lagMult,
                                        PetscScalar temperature,
                                        Vec initialPositions,
                                        PetscScalar NHChainOrder,
                                        PetscInt numIterations,
                                        PetscScalar timestep,
                                        PetscInt numIterationsToSave,
                                        bool saveHamiltonian) /** @todo this might need to be a scalar? */
{
    PetscErrorCode errorStatus = 0;

    if ( (NULL == physics) || (NULL == opt) || (NULL == filter))
    {
        errorStatus = PETSC_ERR_ARG_CORRUPT;
    }

    if (0 == errorStatus)
    {
        /** @todo Make sure all the Vecs are being properly instantiated, this is not correct! */
        this->physics               = physics;
        this->opt                   = opt;
        this->filter                = filter;
        this->lagMult               = lagMult;
        this->temperature           = temperature;
        this->NHChainOrder          = NHChainOrder;
        this->numIterations         = numIterations;
        this->timestep              = timestep;
        this->numIterationsToSave   = numIterationsToSave;
        this->saveHamiltonian       = saveHamiltonian;

        /* Locally set initial values */
        this->numConstraints    = 1; /** @todo This may need to be passed in */
        this->halfTimestep      = timestep/2;

        /* Initialize vectors */
        PetscCall(VecCreate(PETSC_COMM_WORLD, &(this->evenNoseHooverMass)));
        PetscCall(VecSetSizes(this->evenNoseHooverMass, PETSC_DECIDE, NHChainOrder/2));
        PetscCall(VecSetFromOptions(this->evenNoseHooverMass));
        PetscCall(VecDuplicate(this->evenNoseHooverMass, &(this->oddNoseHooverMass)));

        /* */
        PetscCall(VecDuplicate(initialPositions, &(this->newPosition)));
        PetscCall(VecDuplicate(initialPositions, &(this->prevPosition)));
        PetscCall(VecDuplicate(initialPositions, &(this->prevVelocity)));
        PetscCall(VecDuplicate(initialPositions, &(this->sensitivities))); /** @todo Make sure this is correct! */
        PetscCall(VecDuplicate(initialPositions, &(this->constraintSensitivities))); /** @todo Make sure this is correct! */

        /* Set initial values */
        PetscCall(VecSet(this->oddNoseHooverMass,       1.0));
        PetscCall(VecSet(this->evenNoseHooverMass,      1.0));
        PetscCall(VecSet(this->newPosition,             0.0));
        PetscCall(VecSet(this->sensitivities,           1.0)); /** @todo IMPLEMENT*/
        PetscCall(VecSet(this->constraintSensitivities, 1.0)); /** @todo IMPLEMENT*/

        PetscCall(VecCopy(initialPositions, this->prevPosition));
        PetscCall(VecSet(this->prevVelocity, std::sqrt(temperature)));

        PetscInt numPositionParticles;
        PetscCall(VecGetSize(initialPositions, &numPositionParticles));
        this->numParticles = numPositionParticles;

        PetscPrintf(PETSC_COMM_WORLD, "Num Particles: %d\n", this->numParticles);

        if (saveData)
        {
            initializeHDF5();

            LagrangeMultipliers.reserve(numIterations);
            hamiltonians.reserve(numIterations);
            compliance.reserve(numIterations);
            genericData.reserve(numIterations);
            genericData2.reserve(numIterations);
            temperatures.reserve(numIterations);
            iterationTimes.reserve(numIterations);
        }
    }
    return errorStatus;
}

/** @todo FIX THIS */
Hyperoptimization::~Hyperoptimization()
{

}

PetscErrorCode Hyperoptimization::initializeHDF5()
{
    PetscErrorCode errorStatus = 0;

    std::string filename = "hypopt_output";
    std::string fileExtension = ".h5";
    std::string fullPath = "";

    // Check PETSc input for a work directory
    // char      filenameChar[PETSC_MAX_PATH_LEN];
    // PetscBool flg = PETSC_FALSE;
    // PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);

    // If input, change path of the file in filename
    // if (flg) {
    //     fullPath.append(filenameChar);
    //     fullPath.append("/");
    // }
    // else
    // {
    //     fullPath.append("output/");
    // }
    fullPath.append(filename);

    // std::filesystem::path filePath{fullPath + fileExtension};
    // int i = 0;
    std::string appendor = "";
    // while (std::filesystem::exists(filePath))
    // {
    //     i++;
    //     appendor = std::string(i);
    //     filePath = std::filesystem::path(fullPath + appendor + fileExtension);
    // }

    fullPath.append(appendor + fileExtension);

    this->saveFilePath = fullPath;

    PetscPrintf(PETSC_COMM_WORLD, "Trying to save to file: %s\n", fullPath.c_str());

    PetscViewer saveFileHDF5;

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, fullPath.c_str(), FILE_MODE_WRITE, &(saveFileHDF5)));

    // std::time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    // std::string timeString = std::asctime(currentTime);
    // PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/", "create_date", PETSC_STRING, &()));

    PetscCall(PetscViewerHDF5WriteGroup(saveFileHDF5, "/Setting"));
    PetscCall(PetscViewerHDF5WriteGroup(saveFileHDF5, this->stateGroup.c_str()));

    /* Write all settings */

    PetscInt nelx = opt->nxyz[0]-1;
    PetscInt nely = opt->nxyz[1]-1;
    PetscInt nelz = opt->nxyz[2]-1;

    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "volfrac",    PETSC_SCALAR,   &(opt->volfrac)));
    // PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "weights",    PETSC_SCALAR, &(opt->volfrac)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "dt",         PETSC_SCALAR,   &(this->timestep)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "T",          PETSC_SCALAR,   &(this->temperature)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "nelx",       PETSC_INT,      &(nelx)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "nely",       PETSC_INT,      &(nely)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "nelz",       PETSC_INT,      &(nelz)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "penal",      PETSC_SCALAR,   &(opt->penal)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "rmin",       PETSC_SCALAR,   &(opt->rmin)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "stepN",      PETSC_INT,      &(opt->maxItr)));
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "sampleN",    PETSC_INT,      &(opt->maxItr))); /** @todo IMPLEMENT */
    PetscCall(PetscViewerHDF5WriteAttribute(saveFileHDF5, "/Setting", "ch_Order",   PETSC_SCALAR,   &(this->NHChainOrder)));

    PetscCall(PetscViewerDestroy(&(saveFileHDF5)));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::saveIteration(PetscInt iteration, Vec positions)
{
    PetscErrorCode errorStatus = 0;

    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, this->saveFilePath.c_str(), FILE_MODE_APPEND, &(saveFileHDF5)));

    // Set positions name
    std::string stateName = "iteration";
    stateName.append(std::to_string(iteration));

    PetscCall(PetscObjectSetName((PetscObject)positions, stateName.c_str()));
    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, stateGroup.c_str()));
    PetscCall(VecView(positions, saveFileHDF5));
    PetscCall(PetscViewerHDF5PopGroup(saveFileHDF5));

    PetscCall(PetscViewerDestroy(&saveFileHDF5));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::saveFinalValues()
{
    PetscErrorCode errorStatus = 0;

    /* Save Vectors */
    PetscViewer saveFileHDF5;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, this->saveFilePath.c_str(), FILE_MODE_APPEND, &(saveFileHDF5)));
    PetscCall(PetscViewerHDF5PushGroup(saveFileHDF5, dataGroup.c_str()));
    if (saveHamiltonian)
    {
        FileManager::HDF5SaveStdVector(saveFileHDF5, this->hamiltonians, "Hamiltonian");
    }
    PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->temperatures,            "Temperature"));
    PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->LagrangeMultipliers,   "Lambda"));
    PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->compliance,              "Compliance"));
    PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->genericData,             "Volume Fraction"));
    PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->genericData2,            "Max Position"));
    PetscCall(FileManager::HDF5SaveStdVector(saveFileHDF5, this->iterationTimes,           "Iteration Compute Time"));

    PetscCall(PetscObjectSetName((PetscObject)(this->prevVelocity), "Final Velocity"));
    PetscCall(VecView(this->prevVelocity, saveFileHDF5));

    PetscCall(PetscViewerHDF5PopGroup(saveFileHDF5));
    

    PetscCall(PetscViewerDestroy(&saveFileHDF5));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateFirstNoseHooverAcceleration(Vec allVelocities, Vec *accelerations)
{
    PetscErrorCode errorStatus = 0;

    Vec tempVelocities;
    PetscScalar result = 0;
    PetscScalar firstNoseHooverMass;

    /* Get mass */
    PetscCall(PetscExtensions::VecGetOffProcessIndex(this->oddNoseHooverMass, 0, &firstNoseHooverMass));

    /* Square the velocities */
    PetscCall(VecDuplicate(allVelocities, &(tempVelocities)));
    PetscCall(VecCopy(allVelocities, tempVelocities));
    PetscCall(VecPointwiseMult(tempVelocities, tempVelocities, tempVelocities));

    /* Sum over the squared elements */
    PetscCall(VecSum(tempVelocities, &result));

    PetscScalar vsquared = result;
    result = (result - (this->numParticles - this->numConstraints) * this->temperature) / firstNoseHooverMass;

    /**
     * @todo This should be a Vec Copy?
    */
    PetscCall(VecSetValue(*accelerations, 0, result, INSERT_VALUES));

    PetscCall(VecAssemblyBegin(*accelerations));
    PetscCall(VecAssemblyEnd(*accelerations));

    PetscCall(VecDestroy(&tempVelocities));

    return errorStatus;

}

PetscErrorCode Hyperoptimization::calculateRemainingNoseHooverAccelerations(Vec noseHooverVelocities, bool evenOutput, Vec *result)
{
    PetscErrorCode errorStatus = 0;

    Vec decrementedVelocities;
    Vec desiredMasses;
    Vec decrementedMasses;
    Vec tempResults;

    if (evenOutput)
    {
        /* In even case we can use all available Nose Hoover particles */
        PetscCall(VecDuplicate(evenNoseHooverMass, &desiredMasses));
        PetscCall(VecCopy(evenNoseHooverMass, desiredMasses));

        PetscCall(VecDuplicate(oddNoseHooverMass, &decrementedMasses));
        PetscCall(VecCopy(oddNoseHooverMass, decrementedMasses));

        PetscCall(VecDuplicate(noseHooverVelocities, &decrementedVelocities));
        PetscCall(VecCopy(noseHooverVelocities, decrementedVelocities));

        PetscCall(VecDuplicate(desiredMasses, &tempResults));
    }

    else
    {
        PetscInt numReducedParticles = this->NHChainOrder/2 - 1;
        /* In the odd case, we need to omit the first odd mass and last even mass */
        PetscCall(VecCreate(PETSC_COMM_WORLD, &(desiredMasses)));
        PetscCall(VecSetSizes(desiredMasses, PETSC_DETERMINE, numReducedParticles));
        PetscCall(VecSetFromOptions(desiredMasses));

        PetscCall(VecDuplicate(desiredMasses, &decrementedMasses));
        PetscCall(VecDuplicate(desiredMasses, &decrementedVelocities));
        PetscCall(VecDuplicate(desiredMasses, &tempResults));

        Vec sequentialOddNoseHooverMass;
        Vec sequentialEvenNoseHooverMass;
        Vec sequentialNoseHooverVelocities;
        Vec sequentialDesiredMasses;
        Vec sequentialDecrementedMasses;
        Vec sequentialDecrementedVelocities;

        PetscCall(PetscExtensions::VecSequentialFromParallel(this->oddNoseHooverMass,    &sequentialOddNoseHooverMass));
        PetscCall(PetscExtensions::VecSequentialFromParallel(this->evenNoseHooverMass,   &sequentialEvenNoseHooverMass));
        PetscCall(PetscExtensions::VecSequentialFromParallel(noseHooverVelocities,       &sequentialNoseHooverVelocities));
        PetscCall(PetscExtensions::VecSequentialFromParallel(desiredMasses,              &sequentialDesiredMasses));
        PetscCall(PetscExtensions::VecSequentialFromParallel(decrementedMasses,          &sequentialDecrementedMasses));
        PetscCall(PetscExtensions::VecSequentialFromParallel(decrementedVelocities,      &sequentialDecrementedVelocities));

        /* Get the reduced values for each vector */
        const PetscScalar *evenNoseHooverMassArray;
        const PetscScalar *oddNoseHooverMassArray;
        const PetscScalar *noseHooverVelocitiesArray;

        PetscScalar *desiredMassesArray;
        PetscScalar *decrementedMassesArray;
        PetscScalar *decrementedVelocitiesArray;

        PetscCall(VecGetArrayRead(sequentialOddNoseHooverMass,  &oddNoseHooverMassArray));
        PetscCall(VecGetArrayRead(sequentialEvenNoseHooverMass, &evenNoseHooverMassArray));
        PetscCall(VecGetArrayRead(sequentialNoseHooverVelocities,     &noseHooverVelocitiesArray));

        PetscCall(VecGetArray(sequentialDesiredMasses,          &desiredMassesArray));
        PetscCall(VecGetArray(sequentialDecrementedMasses,      &decrementedMassesArray));
        PetscCall(VecGetArray(sequentialDecrementedVelocities,  &decrementedVelocitiesArray));

        for (PetscInt i = 0; i < numReducedParticles; i++)
        {
            desiredMassesArray[i]           = oddNoseHooverMassArray[i+1];
            decrementedMassesArray[i]       = evenNoseHooverMassArray[i];
            decrementedVelocitiesArray[i]   = noseHooverVelocitiesArray[i];
        }

        /* Restore all values */
        PetscCall(VecRestoreArrayRead(sequentialOddNoseHooverMass,  &oddNoseHooverMassArray));
        PetscCall(VecRestoreArrayRead(sequentialEvenNoseHooverMass, &evenNoseHooverMassArray));
        PetscCall(VecRestoreArrayRead(sequentialNoseHooverVelocities,     &noseHooverVelocitiesArray));

        PetscCall(VecRestoreArray(sequentialDesiredMasses,            &desiredMassesArray));
        PetscCall(VecRestoreArray(sequentialDecrementedMasses,        &decrementedMassesArray));
        PetscCall(VecRestoreArray(sequentialDecrementedVelocities,    &decrementedVelocitiesArray));

        PetscCall(PetscExtensions::VecParallelFromSequential(sequentialOddNoseHooverMass,        &(this->oddNoseHooverMass)));
        PetscCall(PetscExtensions::VecParallelFromSequential(sequentialEvenNoseHooverMass,       &(this->evenNoseHooverMass)));
        PetscCall(PetscExtensions::VecParallelFromSequential(sequentialNoseHooverVelocities,     &noseHooverVelocities));
        PetscCall(PetscExtensions::VecParallelFromSequential(sequentialDesiredMasses,            &desiredMasses));
        PetscCall(PetscExtensions::VecParallelFromSequential(sequentialDecrementedMasses,        &decrementedMasses));
        PetscCall(PetscExtensions::VecParallelFromSequential(sequentialDecrementedVelocities,    &decrementedVelocities));

        PetscCall(VecDestroy(&sequentialOddNoseHooverMass));
        PetscCall(VecDestroy(&sequentialEvenNoseHooverMass));
        PetscCall(VecDestroy(&sequentialNoseHooverVelocities));
        PetscCall(VecDestroy(&sequentialDesiredMasses));
        PetscCall(VecDestroy(&sequentialDecrementedMasses));
        PetscCall(VecDestroy(&sequentialDecrementedVelocities));
    }

    /* Calculate the accelerations */
    PetscCall(VecPointwiseMult(decrementedVelocities, decrementedVelocities, decrementedVelocities));
    PetscCall(VecPointwiseMult(tempResults, decrementedMasses, decrementedVelocities));
    PetscCall(VecShift(tempResults, -this->temperature));
    PetscCall(VecPointwiseDivide(tempResults, tempResults, desiredMasses));

    if (evenOutput)
    {
        /* Even output can be directly coppied to the results */
        PetscCall(VecCopy(tempResults, *result));
    }

    else
    {
        Vec resultsSequential;
        Vec tempResultsSequential;

        PetscCall(PetscExtensions::VecSequentialFromParallel(*result, &resultsSequential));
        PetscCall(PetscExtensions::VecSequentialFromParallel(tempResults, &tempResultsSequential));

        /* Skip the first particle for odd results */
        PetscScalar* resultsArray;
        const PetscScalar* tempResultsArray;
        PetscInt numReducedParticles = this->NHChainOrder/2 - 1;

        PetscCall(VecGetArray(resultsSequential, &resultsArray));
        PetscCall(VecGetArrayRead(tempResultsSequential, &tempResultsArray));

        for (PetscInt i = 0; i < numReducedParticles; i++)
        {
            resultsArray[i+1] = tempResultsArray[i];
        }

        PetscCall(VecRestoreArray(resultsSequential, &resultsArray));
        PetscCall(VecRestoreArrayRead(tempResultsSequential, &tempResultsArray));

        PetscCall(PetscExtensions::VecParallelFromSequential(resultsSequential, result));
        PetscCall(PetscExtensions::VecParallelFromSequential(tempResultsSequential, &tempResults));

        PetscCall(VecCopy(*result, *result));

        PetscCall(VecDestroy(&resultsSequential));
        PetscCall(VecDestroy(&tempResultsSequential));

    }

    PetscCall(VecDestroy(&decrementedVelocities));
    PetscCall(VecDestroy(&desiredMasses));
    PetscCall(VecDestroy(&decrementedMasses));
    PetscCall(VecDestroy(&tempResults));

    return errorStatus;
}


PetscErrorCode Hyperoptimization::calculatePositionIncrement(Vec previousPosition, Vec previousVelocity, PetscScalar timeStepIn, Vec *result)
{
    PetscErrorCode errorStatus = 0;

    // Create a temporary vector for inner-results
    Vec tempVector;
    PetscCall(VecDuplicate(previousVelocity, &(tempVector)));
    PetscCall(VecCopy(previousVelocity, tempVector));

    /* dt * v[t-1] */
    PetscCall(VecScale(tempVector, timeStepIn));

    /* x + (dt * v[t-1]) */
    PetscCall(VecAYPX(tempVector, 1, previousPosition));

    /* Copy return value */
    PetscCall(VecCopy(tempVector, *result));

    PetscCall(VecDestroy(&tempVector));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateVelocityIncrement(Vec velocityOne, Vec velocityTwo, Vec acceleration, PetscScalar timeStepIn, Vec *result)
{
    PetscErrorCode errorStatus = 0;

    Vec leftSide;
    Vec rightSide;

    /* Setup left and right vectors */
    PetscCall(VecDuplicate(velocityTwo, &leftSide));
    PetscCall(VecDuplicate(velocityTwo, &rightSide));
    PetscCall(VecCopy(velocityTwo, leftSide));
    PetscCall(VecCopy(velocityTwo, rightSide));

    /* v_1 * exp(- dt * v_2) */
    PetscCall(VecScale(leftSide, - timeStepIn));
    PetscCall(VecExp(leftSide));
    PetscCall(VecPointwiseMult(leftSide, velocityOne, leftSide));

    /* dt * a * exp(- dt/2 * v_2) */
    PetscCall(VecScale(rightSide, - timeStepIn / 2));
    PetscCall(VecExp(rightSide));
    PetscCall(VecPointwiseMult(rightSide, acceleration, rightSide));
    PetscCall(VecScale(rightSide, timeStepIn));

    /* Add left and right and copy to result */
    PetscCall(VecAYPX(leftSide, 1, rightSide));
    PetscCall(VecCopy(leftSide, *result));

    PetscCall(VecDestroy(&leftSide));
    PetscCall(VecDestroy(&rightSide));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateVelocityIncrement(Vec velocityOne, PetscScalar velocityTwo, Vec acceleration, PetscScalar timeStepIn, Vec *result)
{
    PetscErrorCode errorStatus = 0;

    Vec leftSide;
    Vec rightSide;

    /* Setup left and right vectors */
    PetscCall(VecDuplicate(velocityOne, &leftSide));
    PetscCall(VecDuplicate(velocityOne, &rightSide));
    PetscCall(VecCopy(velocityOne, leftSide));
    PetscCall(VecCopy(acceleration, rightSide));

    /* Calculate increment */
    PetscCall(VecScale(leftSide, std::exp( -timeStepIn * velocityTwo )));
    PetscCall(VecScale(rightSide, timeStepIn * std::exp( -timeStepIn * velocityTwo / 2 )));

    /* Add left and right and copy to result */
    PetscCall(VecAYPX(leftSide, 1, rightSide));
    PetscCall(VecCopy(leftSide, *result));

    PetscCall(VecDestroy(&leftSide));
    PetscCall(VecDestroy(&rightSide));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::truncatePositions(Vec *positions)
{
    PetscErrorCode errorStatus = 0;

    PetscInt numPositions;
    PetscCall(VecGetSize(*positions, &numPositions));

    PetscScalar *positionsArray;
    PetscCall(VecGetArray(*positions, &positionsArray));

    PetscInt localSize;
    PetscCall(VecGetLocalSize(*positions, &localSize));

    for (PetscInt i = 0; i < localSize; i++)
    {
        if (positionsArray[i] > 1)
        {
            positionsArray[i] = 1;
        }
        if (positionsArray[i] < 0)
        {
            positionsArray[i] = 0;
        }
    }
    PetscCall(VecRestoreArray(*positions, &positionsArray));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::assembleNewPositions(PetscScalar firstNoseHooverVelocity, PetscScalar *LagrangeMultiplier)
{
    PetscErrorCode errorStatus = 0;

    /* Setup */
    Vec rightSide;
    PetscCall(VecDuplicate(this->newPosition, &(rightSide)));
    PetscCall(VecCopy(this->constraintSensitivities, rightSide));

    /* Lagrange Multipliers require properly constrained design variables */
    truncatePositions(&(this->newPosition));

    // Fix all passive elements
    // opt->SetVariables(this->newPosition, opt->xPassive);

    /* Lagrange multiplier */
    PetscScalar scaleFactor = - this->timestep * this->halfTimestep * std::exp(-this->halfTimestep * firstNoseHooverVelocity);
    PetscCall(VecScale(rightSide, scaleFactor));

    PetscScalar LagrangeMult;
    this->lagMult.computeLagrangeMultiplier(this->newPosition, rightSide, this->numParticles, &LagrangeMult);

    PetscCall(VecScale(rightSide, LagrangeMult));
    PetscCall(VecAYPX(this->newPosition, 1, rightSide));

    *LagrangeMultiplier = LagrangeMult;

    truncatePositions(&(this->newPosition));

    // Fix all passive elements
    // opt->SetVariables(this->newPosition, opt->xPassive);

    PetscCall(VecDestroy(&rightSide));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateSensitvities(Vec positions)
{
    PetscErrorCode errorStatus = 0;
    // Positions used in constraint calculations
    errorStatus = this->filter->FilterProject(positions, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);
    CHKERRQ(errorStatus);

    // Compute sensitivities
    errorStatus = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), 
                                                                    &(opt->gx[0]),
                                                                    opt->dfdx,
                                                                    opt->dgdx[0],
                                                                    opt->xTilde,
                                                                    // opt->xPhys,
                                                                    opt->Emin,
                                                                    opt->Emax,
                                                                    opt->penal,
                                                                    opt->volfrac);//,
                                                                    // data);
    CHKERRQ(errorStatus);

    // Filter sensitivities (chainrule)
    errorStatus = filter->Gradients(opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, opt->beta,
                            opt->eta);
    CHKERRQ(errorStatus);

    PetscCall(VecCopy(opt->dfdx, this->sensitivities));
    PetscCall(VecScale(this->sensitivities, -1));
    PetscCall(VecCopy(opt->dgdx[0], this->constraintSensitivities));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateTemperature(Vec velocities, PetscScalar *systemTemperature)
{
    PetscErrorCode errorStatus = 0;

    Vec tempVector;
    PetscCall(VecDuplicate(velocities, &tempVector));
    PetscCall(VecCopy(velocities, tempVector));

    PetscCall(VecPointwiseMult(tempVector, tempVector, tempVector));
    PetscCall(VecMean(tempVector, systemTemperature));

    PetscCall(VecDestroy(&tempVector));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateHamiltonian(Vec velocities, Vec positions, PetscScalar *hamiltonian)
{
    PetscErrorCode errorStatus = 0;

    PetscScalar velocityDotProduct;

    truncatePositions(&positions);

    errorStatus = this->filter->FilterProject(positions, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);

    errorStatus = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), 
                                                                    &(opt->gx[0]),
                                                                    opt->dfdx,
                                                                    opt->dgdx[0],
                                                                    opt->xTilde,
                                                                    // opt->xPhys,
                                                                    opt->Emin,
                                                                    opt->Emax,
                                                                    opt->penal,
                                                                    opt->volfrac);//,

    PetscScalar currentCompliance = opt->fx;

    PetscCall(VecDot(velocities, velocities, &velocityDotProduct));

    *hamiltonian = currentCompliance + velocityDotProduct / 2;

    this->compliance.push_back(currentCompliance);

    return errorStatus;
}

PetscErrorCode Hyperoptimization::runDesignLoop()
{
    PetscErrorCode errorStatus = 0;

    PetscPrintf(PETSC_COMM_WORLD, "Initializing...");

    PetscInt numOddNHIndices = this->NHChainOrder/2;
    PetscInt numEvenNHIndices = this->NHChainOrder/2;
    PetscInt straightNHIndices[numOddNHIndices];

    /* NHChainOrder is assumed even, confirm with Hazhir if this is true */
    for (PetscInt i = 0; i < numEvenNHIndices; i++)
    {
        straightNHIndices[i] = i;
    }

    Vec tempPosition;
    // Vec prevPosition;
    Vec newVelocity;

    Vec tempOddNoseHooverVelocity;
    Vec tempEvenNoseHooverPosition;
    Vec tempOddNoseHooverAccelerations;
    Vec tempEvenNoseHooverAccelerations;

    Vec prevEvenNoseHooverPosition;
    Vec prevEvenNoseHooverVelocity;
    Vec prevOddNoseHooverPosition;
    Vec prevOddNoseHooverVelocity;

    Vec newOddNoseHooverPosition;
    Vec newEvenNoseHooverPosition;
    Vec newEvenNoseHooverVelocity;
    Vec newOddNoseHooverVelocity;
    // Vec extendedOddNoseHooverVelocity;

    /* Setup all vectors */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(tempOddNoseHooverVelocity)));
    PetscCall(VecSetSizes(tempOddNoseHooverVelocity, PETSC_DETERMINE, numOddNHIndices));
    PetscCall(VecSetFromOptions(tempOddNoseHooverVelocity));

    PetscCall(VecDuplicate(this->newPosition, &(tempPosition)));
    PetscCall(VecDuplicate(this->newPosition, &(newVelocity)));

    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(tempOddNoseHooverVelocity)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(tempEvenNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(tempOddNoseHooverAccelerations)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(tempEvenNoseHooverAccelerations)));

    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(prevEvenNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(prevEvenNoseHooverVelocity)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(prevOddNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(prevOddNoseHooverVelocity)));

    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newOddNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newEvenNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newEvenNoseHooverVelocity)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newOddNoseHooverVelocity)));
    // PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &extendedOddNoseHooverVelocity));

    /* Initial Values */
    PetscCall(VecSet(prevEvenNoseHooverPosition,    0));
    PetscCall(VecSet(prevEvenNoseHooverVelocity,    0));
    PetscCall(VecSet(prevOddNoseHooverPosition,     0));
    PetscCall(VecSet(prevOddNoseHooverVelocity,     0));

    PetscPrintf(PETSC_COMM_WORLD, "...initialized!\n");

    for (PetscInt iteration = 0; iteration < this->numIterations; iteration++)
    {
        double t1 = MPI_Wtime();

        /* Calculate half-timestep positions */
        calculatePositionIncrement(prevPosition, this->prevVelocity, this->halfTimestep, &tempPosition);
        calculatePositionIncrement(prevEvenNoseHooverPosition, prevEvenNoseHooverVelocity, this->halfTimestep, &tempEvenNoseHooverPosition);

        /* Calcualte previous Nose Hoover accelerations */
        calculateFirstNoseHooverAcceleration(this->prevVelocity, &tempOddNoseHooverAccelerations);
        calculateRemainingNoseHooverAccelerations(prevEvenNoseHooverVelocity, false, &tempOddNoseHooverAccelerations);

        /* Calculate half-timestep velocities */
        calculateVelocityIncrement(prevOddNoseHooverVelocity, prevEvenNoseHooverVelocity, tempOddNoseHooverAccelerations, this->halfTimestep, &tempOddNoseHooverVelocity);

        /* Get objective and constraint sensitivities */
        calculateSensitvities(tempPosition);

        /* Get first nose hoover velocity */
        PetscScalar firstNoseHooverVelocity;
        const PetscScalar *tempOddNoseHooverVelocityArray;

        PetscCall(PetscExtensions::VecGetOffProcessIndex(tempOddNoseHooverVelocity, 0, &firstNoseHooverVelocity));

        /* Calculate full increment values */
        calculateVelocityIncrement(this->prevVelocity, firstNoseHooverVelocity, this->sensitivities, this->timestep, &newVelocity);
        calculatePositionIncrement(prevOddNoseHooverPosition, tempOddNoseHooverVelocity, this->timestep, &newOddNoseHooverPosition);

        /* Nose Hoover acceleration at half-timestep */
        calculateRemainingNoseHooverAccelerations(tempOddNoseHooverVelocity, true, &tempEvenNoseHooverAccelerations);

        /* Add 0 to the end of the odd nose hoover velocities */
        Vec extendedOddNoseHooverVelocity;
        PetscExtensions::VecLeftShift(tempOddNoseHooverVelocity, &extendedOddNoseHooverVelocity);

        /* Continue calculating full-timestep values */
        calculateVelocityIncrement(prevEvenNoseHooverVelocity, extendedOddNoseHooverVelocity, tempEvenNoseHooverAccelerations, this->timestep, &(newEvenNoseHooverVelocity));
        calculatePositionIncrement(tempPosition, newVelocity, this->halfTimestep, &(this->newPosition));
        calculatePositionIncrement(tempEvenNoseHooverPosition, newEvenNoseHooverVelocity, this->halfTimestep, &(newEvenNoseHooverPosition));

        /* Calculate the new positions */
        PetscScalar lagMultiplier;
        assembleNewPositions(firstNoseHooverVelocity, &lagMultiplier);

        /* Calcualte the new velocities as v_new = (x_new -  x_old) / timestep */
        PetscCall(VecWAXPY(newVelocity, -1, prevPosition, newPosition));
        PetscCall(VecScale(newVelocity, 1/this->timestep));

        /* Calculate full-timestep Nose Hoover accelerations */
        calculateFirstNoseHooverAcceleration(newVelocity, &tempOddNoseHooverAccelerations);
        calculateRemainingNoseHooverAccelerations(newEvenNoseHooverVelocity, false, &tempOddNoseHooverAccelerations);

        /* Calculate final odd Nose Hoover Velocity */
        calculateVelocityIncrement(tempOddNoseHooverVelocity, newEvenNoseHooverVelocity, tempOddNoseHooverAccelerations, this->halfTimestep, &newOddNoseHooverVelocity);

        double t2 = MPI_Wtime();

        if (saveData)
        {
            if ( (this->numIterations - iteration) <= numIterationsToSave)
            {
                saveIteration(iteration, this->newPosition);
            }

            if (temperatureCheck)
            {
                Vec filtered_pos;

                PetscCall(VecDuplicate(this->newPosition, &filtered_pos));
                errorStatus = this->filter->FilterProject(newPosition, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);

                PetscCall(VecCopy(opt->xTilde, filtered_pos));

                PetscScalar temperature;
                calculateTemperature(newVelocity, &temperature);

                PetscReal minPos;
                PetscReal maxPos;
                PetscScalar meanPos;

                PetscCall(VecMax(filtered_pos, NULL, &maxPos));
                PetscCall(VecMin(filtered_pos, NULL, &minPos));
                PetscCall(VecMean(filtered_pos, &meanPos));

                PetscReal minVel;
                PetscReal maxVel;
                PetscScalar meanVel;

                PetscCall(VecMax(newVelocity, NULL, &maxVel));
                PetscCall(VecMin(newVelocity, NULL, &minVel));
                PetscCall(VecMean(newVelocity, &meanVel));

                if (saveHamiltonian)
                {
                    PetscScalar hamiltonian;
                    calculateHamiltonian(newVelocity, this->newPosition, &hamiltonian);
                    hamiltonians.push_back(hamiltonian);
                    t2 = MPI_Wtime();
                }

                temperatures.push_back(temperature);
                LagrangeMultipliers.push_back(lagMultiplier);
                genericData.push_back(meanPos);
                genericData2.push_back(maxPos);
                iterationTimes.push_back(t2 - t1);

                // PetscPrintf(PETSC_COMM_WORLD, "iter: %i, Max Pos: %f, Min Pos: %f, Mean Pos: %f, Max Vel: %f, Min Vel: %f, Mean Vel: %f, Temperature: %f, LM: %f, HM: %f\n", iteration, maxPos, minPos, meanPos, maxVel, minVel, meanVel, temperature, lagMultiplier);//, hamiltonian);
                PetscPrintf(PETSC_COMM_WORLD, "iter: %i, Max Pos: %f, Min Pos: %f, Mean Pos: %f, Max Vel: %f, Min Vel: %f, Mean Vel: %f, Temp: %f, LM: %f,\ttime: %f\n", iteration, maxPos, minPos, meanPos, maxVel, minVel, meanVel, temperature, lagMultiplier, t2 - t1);//, hamiltonian);

                PetscCall(VecDestroy(&filtered_pos));
            }
        }

        /* SAVE VARIABLES! */
        PetscCall(VecCopy(newPosition, prevPosition));
        PetscCall(VecCopy(newVelocity, prevVelocity));

        PetscCall(VecCopy(newEvenNoseHooverPosition,    prevEvenNoseHooverPosition));
        PetscCall(VecCopy(newEvenNoseHooverVelocity,    prevEvenNoseHooverVelocity));
        PetscCall(VecCopy(newOddNoseHooverPosition,     prevOddNoseHooverPosition));
        PetscCall(VecCopy(newOddNoseHooverVelocity,     prevOddNoseHooverVelocity));
    }

    if (saveData)
    {
        saveFinalValues();
    }

    /* Cleanup! */
    PetscCall(VecDestroy(&tempPosition));
    PetscCall(VecDestroy(&newVelocity));
    PetscCall(VecDestroy(&tempOddNoseHooverVelocity));
    PetscCall(VecDestroy(&tempEvenNoseHooverPosition));
    PetscCall(VecDestroy(&tempOddNoseHooverAccelerations));
    PetscCall(VecDestroy(&tempEvenNoseHooverAccelerations));
    PetscCall(VecDestroy(&prevEvenNoseHooverPosition));
    PetscCall(VecDestroy(&prevEvenNoseHooverVelocity));
    PetscCall(VecDestroy(&prevOddNoseHooverPosition));
    PetscCall(VecDestroy(&prevOddNoseHooverVelocity));
    PetscCall(VecDestroy(&newOddNoseHooverPosition));
    PetscCall(VecDestroy(&newEvenNoseHooverPosition));
    PetscCall(VecDestroy(&newEvenNoseHooverVelocity));
    PetscCall(VecDestroy(&newOddNoseHooverVelocity));
    // PetscCall(VecDestroy(&extendedOddNoseHooverVelocity));

    return errorStatus;
}
