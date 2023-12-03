/***************************************************************************//**
 * @file Hyperoptimization.cc
 *
 * This file describes the layout of the Hyperoptimization class.
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#include "Hyperoptimization.h"
#include "HypOptException.h"

#include "PetscExtensions.h"
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
 * 
 * @todo Convert ALL wrappers to use pass by reference using &!!!!!! This will ensure that overwritten functions
 * are actually processed properly.
*/
PetscErrorCode Hyperoptimization::init( SensitivitiesWrapper&   sensitivitiesWrapper,
                                        FilterWrapper&          filter,
                                        LagrangeMultiplier&     lagMult,
                                        PetscScalar             temperature,
                                        Vec                     initialPositions,
                                        Vec                     initialVelocities,
                                        PetscInt                NHChainOrder,
                                        PetscInt                numIterations,
                                        PetscScalar             timestep,
                                        FileManager*            fileManager,
                                        std::vector<uint32_t>   iterationSaveRange,
                                        bool                    saveHamiltonian) /** @todo this might need to be a scalar? */
{
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(prevState.evenNoseHooverPosition)));
    PetscCall(VecSetSizes(prevState.evenNoseHooverPosition, PETSC_DETERMINE, NHChainOrder/2));
    PetscCall(VecSetFromOptions(prevState.evenNoseHooverPosition));

    PetscCall(VecDuplicate(prevState.evenNoseHooverPosition, &(prevState.evenNoseHooverVelocity)));
    PetscCall(VecDuplicate(prevState.evenNoseHooverPosition, &(prevState.oddNoseHooverPosition)));
    PetscCall(VecDuplicate(prevState.evenNoseHooverPosition, &(prevState.oddNoseHooverVelocity)));

    /* initial Values */
    PetscCall(VecSet(prevState.evenNoseHooverPosition,    0));
    PetscCall(VecSet(prevState.evenNoseHooverVelocity,    0));
    PetscCall(VecSet(prevState.oddNoseHooverPosition,     0));
    PetscCall(VecSet(prevState.oddNoseHooverVelocity,     0));

    return init(sensitivitiesWrapper,
                lagMult,
                filter,
                fileManager,
                NHChainOrder,
                timestep,
                temperature,
                numIterations,
                iterationSaveRange,
                initialPositions,
                initialVelocities,
                prevState.evenNoseHooverPosition,
                prevState.evenNoseHooverVelocity,
                prevState.oddNoseHooverPosition,
                prevState.oddNoseHooverVelocity,
                saveHamiltonian);
}

PetscErrorCode Hyperoptimization::init( SensitivitiesWrapper&   sensitivitiesWrapper,
                                        LagrangeMultiplier&     lagMult,
                                        FilterWrapper&          filter,
                                        FileManager*            fileManager,
                                        PetscInt                NHChainOrder,
                                        PetscScalar             timestep,
                                        PetscScalar             temperature,
                                        PetscInt                numIterations,
                                        std::vector<uint32_t>   iterationSaveRange,
                                        Vec                     initialPositions,
                                        Vec                     initialVelocities,
                                        Vec                     initialEvenNoseHooverPosition,
                                        Vec                     initialEvenNoseHooverVelocity,
                                        Vec                     initialOddNoseHooverPosition,
                                        Vec                     initialOddNoseHooverVelocity,
                                        bool                    saveHamiltonian) /** @todo this might need to be a scalar? */
{
    PetscErrorCode errorStatus = 0;

    if (NULL == fileManager)
    {
        throw HypOptException("Null pointer error.");
    }

    if (2 != iterationSaveRange.size()          ||
        iterationSaveRange[0] < 0               ||
        iterationSaveRange[1] > numIterations   ||
        iterationSaveRange[0] > iterationSaveRange[1])
    {
        throw HypOptException("Invalid iteration save range. Must be two integers between 0 and the number of iterations, the first lower than the second.");
    }

    /** @todo Make sure all the Vecs are being properly instantiated, this is not correct! */
    this->fileManager           = fileManager;
    this->sensitivitiesWrapper          = sensitivitiesWrapper;
    this->filter                = filter;
    this->lagMult               = lagMult;
    this->temperature           = temperature;
    this->NHChainOrder          = NHChainOrder;
    this->numIterations         = numIterations;
    this->timestep              = timestep;
    this->iterationSaveRange    = iterationSaveRange;
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
    PetscCall(VecDuplicate(initialPositions,                &(this->newPosition)));
    PetscCall(VecDuplicate(initialPositions,                &(this->prevState.position)));
    PetscCall(VecDuplicate(initialVelocities,               &(this->prevState.velocity)));
    PetscCall(VecDuplicate(initialEvenNoseHooverPosition,   &(this->prevState.evenNoseHooverPosition)));
    PetscCall(VecDuplicate(initialEvenNoseHooverVelocity,   &(this->prevState.evenNoseHooverVelocity)));
    PetscCall(VecDuplicate(initialOddNoseHooverPosition,    &(this->prevState.oddNoseHooverPosition)));
    PetscCall(VecDuplicate(initialOddNoseHooverVelocity,    &(this->prevState.oddNoseHooverVelocity)));
    PetscCall(VecDuplicate(initialPositions,                &(this->sensitivities))); /** @todo Make sure this is correct! */
    PetscCall(VecDuplicate(initialPositions,                &(this->constraintSensitivities))); /** @todo Make sure this is correct! */

    /* Set initial values */
    PetscCall(VecSet(this->oddNoseHooverMass,       1.0));
    PetscCall(VecSet(this->evenNoseHooverMass,      1.0));
    PetscCall(VecSet(this->newPosition,             0.0));
    PetscCall(VecSet(this->sensitivities,           1.0)); /** @todo IMPLEMENT*/
    PetscCall(VecSet(this->constraintSensitivities, 1.0)); /** @todo IMPLEMENT*/

    PetscCall(VecCopy(initialPositions,                 this->prevState.position));
    PetscCall(VecCopy(initialVelocities,                this->prevState.velocity));
    PetscCall(VecCopy(initialEvenNoseHooverPosition,    this->prevState.evenNoseHooverPosition));
    PetscCall(VecCopy(initialEvenNoseHooverVelocity,    this->prevState.evenNoseHooverVelocity));
    PetscCall(VecCopy(initialOddNoseHooverPosition,     this->prevState.oddNoseHooverPosition));
    PetscCall(VecCopy(initialOddNoseHooverVelocity,     this->prevState.oddNoseHooverVelocity));

    PetscInt numPositionParticles;
    PetscCall(VecGetSize(initialPositions, &numPositionParticles));
    this->numParticles = numPositionParticles;

    if (saveData)
    {
        LagrangeMultipliers.reserve(numIterations);
        temperatures.reserve(numIterations);
        iterationTimes.reserve(numIterations);

        if (saveHamiltonian)
        {
            hamiltonians.reserve(numIterations);
            compliance.reserve(numIterations);
        }
    }

 /** @todo print out hypopt settings */

    PetscPrintf(PETSC_COMM_WORLD, "#\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Hyperoptimization setup. Settings:\n");
    PetscPrintf(PETSC_COMM_WORLD, "# -targetTemperature: %f\n", temperature);
    PetscPrintf(PETSC_COMM_WORLD, "# -NHChainOrder: %i\n", NHChainOrder);
    PetscPrintf(PETSC_COMM_WORLD, "# -maximumIterations: %i\n", numIterations);
    PetscPrintf(PETSC_COMM_WORLD, "# -timestep: %f\n", timestep);
    PetscPrintf(PETSC_COMM_WORLD, "# -iterationSaveRange: (%i, %i)\n", iterationSaveRange[0], iterationSaveRange[1]);

    return errorStatus;
}

/** @todo FIX THIS */
Hyperoptimization::~Hyperoptimization()
{

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

    /* Lagrange multiplier */
    PetscScalar scaleFactor = - this->timestep * this->halfTimestep * std::exp(-this->halfTimestep * firstNoseHooverVelocity);
    PetscCall(VecScale(rightSide, scaleFactor));

    PetscScalar LagrangeMult;
    this->lagMult.computeLagrangeMultiplier(this->newPosition, rightSide, this->numParticles, &LagrangeMult);

    PetscCall(VecScale(rightSide, LagrangeMult));
    PetscCall(VecAYPX(this->newPosition, 1, rightSide));

    if (nullptr != LagrangeMultiplier)
    {
        *LagrangeMultiplier = LagrangeMult;
    }

    truncatePositions(&(this->newPosition));

    PetscCall(VecDestroy(&rightSide));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateSensitvities(Vec positions)
{
    PetscErrorCode errorStatus = 0;

    Vec filteredPositions;
    PetscCall(VecDuplicate(positions, &filteredPositions));

    // Positions used in constraint calculations
    PetscCall(this->filter.filterDesignVariable(positions, filteredPositions));

    // Compute sensitivities
    PetscCall(sensitivitiesWrapper.computeSensitivities(filteredPositions, sensitivities, constraintSensitivities));

    // Filter sensitivities (Standard Filter)
    PetscCall(filter.filterSensitivities(positions, this->sensitivities, &constraintSensitivities));
    PetscCall(VecScale(this->sensitivities, -1));

    PetscCall(VecDestroy(&filteredPositions));

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
    PetscScalar currentCompliance = 0;
    PetscScalar velocityDotProduct;

    Vec filteredPositions;
    PetscCall(VecDuplicate(positions, &filteredPositions));

    truncatePositions(&positions);
    errorStatus = filter.filterDesignVariable(positions, filteredPositions);
    errorStatus = sensitivitiesWrapper.computeObjectiveFunction(filteredPositions, &currentCompliance);
    PetscCall(VecDot(velocities, velocities, &velocityDotProduct));
    *hamiltonian = currentCompliance + velocityDotProduct / 2;

    this->compliance.push_back(currentCompliance);

    PetscCall(VecDestroy(&filteredPositions));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::runDesignLoop()
{
    PetscErrorCode errorStatus = 0;

    PetscInt numOddNHIndices = this->NHChainOrder/2;
    PetscInt numEvenNHIndices = this->NHChainOrder/2;

    Vec tempPosition;
    // Vec prevState.position;
    Vec newVelocity;

    Vec tempOddNoseHooverVelocity;
    Vec tempEvenNoseHooverPosition;
    Vec tempOddNoseHooverAccelerations;
    Vec tempEvenNoseHooverAccelerations;

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

    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newOddNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newEvenNoseHooverPosition)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newEvenNoseHooverVelocity)));
    PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &(newOddNoseHooverVelocity)));
    // PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &extendedOddNoseHooverVelocity));

    for (PetscInt iteration = 0; iteration < this->numIterations; iteration++)
    {
        double t1 = MPI_Wtime();

        /* Calculate half-timestep positions */
        calculatePositionIncrement(prevState.position, this->prevState.velocity, this->halfTimestep, &tempPosition);
        calculatePositionIncrement(prevState.evenNoseHooverPosition, prevState.evenNoseHooverVelocity, this->halfTimestep, &tempEvenNoseHooverPosition);

        /* Calcualte previous Nose Hoover accelerations */
        calculateFirstNoseHooverAcceleration(this->prevState.velocity, &tempOddNoseHooverAccelerations);
        calculateRemainingNoseHooverAccelerations(prevState.evenNoseHooverVelocity, false, &tempOddNoseHooverAccelerations);

        /* Calculate half-timestep velocities */
        calculateVelocityIncrement(prevState.oddNoseHooverVelocity, prevState.evenNoseHooverVelocity, tempOddNoseHooverAccelerations, this->halfTimestep, &tempOddNoseHooverVelocity);

        /* Get objective and constraint sensitivities */
        calculateSensitvities(tempPosition);

        /* Get first nose hoover velocity */
        PetscScalar firstNoseHooverVelocity;
        const PetscScalar *tempOddNoseHooverVelocityArray;

        PetscCall(PetscExtensions::VecGetOffProcessIndex(tempOddNoseHooverVelocity, 0, &firstNoseHooverVelocity));

        /* Calculate full increment values */
        calculateVelocityIncrement(this->prevState.velocity, firstNoseHooverVelocity, this->sensitivities, this->timestep, &newVelocity);
        calculatePositionIncrement(prevState.oddNoseHooverPosition, tempOddNoseHooverVelocity, this->timestep, &newOddNoseHooverPosition);

        /* Nose Hoover acceleration at half-timestep */
        calculateRemainingNoseHooverAccelerations(tempOddNoseHooverVelocity, true, &tempEvenNoseHooverAccelerations);

        /* Add 0 to the end of the odd nose hoover velocities */
        Vec extendedOddNoseHooverVelocity;
        PetscExtensions::VecLeftShift(tempOddNoseHooverVelocity, &extendedOddNoseHooverVelocity);

        /* Continue calculating full-timestep values */
        calculateVelocityIncrement(prevState.evenNoseHooverVelocity, extendedOddNoseHooverVelocity, tempEvenNoseHooverAccelerations, this->timestep, &(newEvenNoseHooverVelocity));
        calculatePositionIncrement(tempPosition, newVelocity, this->halfTimestep, &(this->newPosition));
        calculatePositionIncrement(tempEvenNoseHooverPosition, newEvenNoseHooverVelocity, this->halfTimestep, &(newEvenNoseHooverPosition));

        /* Calculate the new positions */
        PetscScalar lagMultiplier;
        assembleNewPositions(firstNoseHooverVelocity, &lagMultiplier);

        /* Calcualte the new velocities as v_new = (x_new -  x_old) / timestep */
        PetscCall(VecWAXPY(newVelocity, -1, prevState.position, newPosition));
        PetscCall(VecScale(newVelocity, 1/this->timestep));

        /* Calculate full-timestep Nose Hoover accelerations */
        calculateFirstNoseHooverAcceleration(newVelocity, &tempOddNoseHooverAccelerations);
        calculateRemainingNoseHooverAccelerations(newEvenNoseHooverVelocity, false, &tempOddNoseHooverAccelerations);

        /* Calculate final odd Nose Hoover Velocity */
        calculateVelocityIncrement(tempOddNoseHooverVelocity, newEvenNoseHooverVelocity, tempOddNoseHooverAccelerations, this->halfTimestep, &newOddNoseHooverVelocity);

        double t2 = MPI_Wtime();

        if (saveData)
        {
            if ( (this->numIterations - iteration) <= iterationSaveRange[1] &&
                 (this->numIterations - iteration) >= iterationSaveRange[0]  )
            {
                fileManager->saveIteration(iteration, this->newPosition);
            }

            if (temperatureCheck)
            {
                Vec filtered_pos;

                PetscCall(VecDuplicate(this->newPosition, &filtered_pos));
                errorStatus = filter.filterDesignVariable(newPosition, filtered_pos);

                PetscScalar systemTemperature;
                calculateTemperature(newVelocity, &systemTemperature);

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

                temperatures.push_back(systemTemperature);
                LagrangeMultipliers.push_back(lagMultiplier);
                iterationTimes.push_back(t2 - t1);

                // PetscPrintf(PETSC_COMM_WORLD, "iter: %i, Max Pos: %f, Min Pos: %f, Mean Pos: %f, Max Vel: %f, Min Vel: %f, Mean Vel: %f, Temperature: %f, LM: %f, HM: %f\n", iteration, maxPos, minPos, meanPos, maxVel, minVel, meanVel, temperature, lagMultiplier);//, hamiltonian);
                PetscPrintf(PETSC_COMM_WORLD, "iter: %i, Max Pos: %e, Min Pos: %e, Mean Pos: %e, Max Vel: %e, Min Vel: %e, Mean Vel: %e, Temp: %e, LM: %e,\ttime: %f\n", iteration, maxPos, minPos, meanPos, maxVel, minVel, meanVel, systemTemperature, lagMultiplier, t2 - t1);//, hamiltonian);

                PetscCall(VecDestroy(&filtered_pos));
            }
        }

        /* SAVE VARIABLES! */
        PetscCall(VecCopy(newPosition, prevState.position));
        PetscCall(VecCopy(newVelocity, prevState.velocity));

        PetscCall(VecCopy(newEvenNoseHooverPosition,    prevState.evenNoseHooverPosition));
        PetscCall(VecCopy(newEvenNoseHooverVelocity,    prevState.evenNoseHooverVelocity));
        PetscCall(VecCopy(newOddNoseHooverPosition,     prevState.oddNoseHooverPosition));
        PetscCall(VecCopy(newOddNoseHooverVelocity,     prevState.oddNoseHooverVelocity));
    }

    doneSolving = true;

    /* Cleanup! */
    PetscCall(VecDestroy(&tempPosition));
    PetscCall(VecDestroy(&newVelocity));
    PetscCall(VecDestroy(&tempOddNoseHooverVelocity));
    PetscCall(VecDestroy(&tempEvenNoseHooverPosition));
    PetscCall(VecDestroy(&tempOddNoseHooverAccelerations));
    PetscCall(VecDestroy(&tempEvenNoseHooverAccelerations));
    PetscCall(VecDestroy(&newOddNoseHooverPosition));
    PetscCall(VecDestroy(&newEvenNoseHooverPosition));
    PetscCall(VecDestroy(&newEvenNoseHooverVelocity));
    PetscCall(VecDestroy(&newOddNoseHooverVelocity));
    // PetscCall(VecDestroy(&extendedOddNoseHooverVelocity));

    return errorStatus;
}
