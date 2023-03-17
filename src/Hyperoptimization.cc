/**
 * @file Hyperoptimization.h
 * 
 * This file describes the layout of the Hyperoptimization class.
 * 
 * @todo LICENSE
**/

#include "Hyperoptimization.h"

#include <cmath>


/**
 * @todo ALL CALLS TO VecSetValues MUST be followed up by "VecAssemblyBegin()" and "VecAssemblyEnd()"!!!!
 *       GO THROUGH THE WHOLE CODE AND UPDATE THIS.
 * 
 * @todo ALL CALLS TO VecGetArrayRead MUST be followed up by "VecRestoreArrayRead()"!!!!!!!!!
*/

// Hyperoptimization::Hyperoptimization(   LinearElasticity* physics,
//                                         TopOpt* opt,
//                                         Filter* filter,
//                                         DataObj data,
//                                         LagrangianMultiplier lagMult,
//                                         PetscScalar temperature,
//                                         Vec initialPositions,
//                                         Vec initialVelocities,
//                                         PetscScalar NHChainOrder,
//                                         PetscInt numIterations,
//                                         PetscFloat timestep)
// {}

PetscErrorCode Hyperoptimization::init( LinearElasticity* physics,
                                        TopOpt* opt,
                                        Filter* filter,
                                        // DataObj data,
                                        LagrangianMultiplier lagMult,
                                        PetscScalar temperature,
                                        Vec initialPositions,
                                        PetscScalar NHChainOrder,
                                        PetscInt numIterations,
                                        PetscFloat timestep) /** @todo this might need to be a scalar? */
{
    PetscErrorCode errorStatus;

    if (nullptr == physics)
    {
        // Throw exception (or petsc equivalent??)
    }

    if (nullptr == opt)
    {
        // Throw exception (or petsc equivalent??)
    }

    /** @todo Make sure all the Vecs are being properly instantiated, this is not correct! */
    this->physics           = physics;
    this->opt               = opt;
    this->filter            = filter;
    // this->data              = data;
    this->lagMult           = lagMult;
    this->temperature       = temperature;
    this->NHChainOrder      = NHChainOrder;
    this->numIterations     = numIterations;
    this->timestep          = timestep;


    /* Locally set initial values */
    this->numConstraints    = 1; /** @todo This may need to be passed in */
    this->halfTimestep      = timestep/2;

    /* Initialize vectors */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(this->evenNoseHooverMass)));
    PetscCall(VecSetSizes(this->evenNoseHooverMass, PETSC_DECIDE, NHChainOrder/2));
    PetscCall(VecSetFromOptions(this->evenNoseHooverMass));
    PetscCall(VecDuplicate(this->evenNoseHooverMass, &(this->oddNoseHooverMass)));

    PetscCall(VecDuplicate(initialPositions, &(this->newPosition)));
    PetscCall(VecDuplicate(initialPositions, &(this->prevPosition)));
    PetscCall(VecDuplicate(initialPositions, &(this->prevVelocity)));
    PetscCall(VecDuplicate(initialPositions, &(this->sensitivities))); /** @todo Make sure this is correct! */
    PetscCall(VecDuplicate(initialPositions, &(this->constraintSensitivities))); /** @todo Make sure this is correct! */

    /** @todo Setup structure to save all iterations */
    // this->positions.resize(numIterations + 1);

    // for (Vec iteration : positions)
    // {
    //     PetscCall(VecDuplicate(initialPositions, &iteration));
    // }

    /* Set initial values */
    PetscCall(VecSet(this->oddNoseHooverMass,       1.0));
    PetscCall(VecSet(this->evenNoseHooverMass,      1.0));
    PetscCall(VecSet(this->newPosition,             0.0));
    PetscCall(VecSet(this->sensitivities,           1.0)); /** @todo IMPLEMENT*/
    PetscCall(VecSet(this->constraintSensitivities, 1.0)); /** @todo IMPLEMENT*/

    PetscCall(VecCopy(initialPositions, this->prevPosition));

    // Temperature = mean (velocities^2) -> initial velocities = sqrt(temperature)
    PetscCall(VecSet(this->prevVelocity, std::sqrt(temperature)));

    PetscInt numPositionParticles;
    PetscCall(VecGetSize(initialPositions, &numPositionParticles));
    this->numParticles = numPositionParticles;

    PetscPrintf(PETSC_COMM_WORLD, "Num Particles: %d\n", this->numParticles);

    return errorStatus;
}

/** @todo FIX THIS */
Hyperoptimization::~Hyperoptimization()
{
    // if (NULL != initialPositions)
    // {
    //     VecDestroy(&initialPositions);
    // }

    // if (NULL != initialVelocities)
    // {
    //     VecDestroy(&initialVelocities);
    // }

    // if (NULL != passiveElements)
    // {
    //     VecDestroy(&passiveElements);
    // }

}


PetscErrorCode Hyperoptimization::calculateFirstNoseHooverAcceleration(Vec allVelocities, Vec *accelerations)
{
    PetscErrorCode errorStatus;

    Vec tempVelocities;
    PetscScalar result = 0;
    PetscScalar firstNoseHooverMass;

    /* Get mass */
    PetscInt index = 0;
    const PetscScalar *noseHooverMassArray;

    PetscCall(VecGetArrayRead(this->oddNoseHooverMass, &noseHooverMassArray));
    firstNoseHooverMass = noseHooverMassArray[0];
    PetscCall(VecRestoreArrayRead(this->oddNoseHooverMass, &noseHooverMassArray));

    /* Square the velocities */
    PetscCall(VecDuplicate(allVelocities, &(tempVelocities)));
    PetscCall(VecCopy(allVelocities, tempVelocities));
    PetscCall(VecPointwiseMult(tempVelocities, tempVelocities, tempVelocities));

    /* Sum over the squared elements */
    PetscCall(VecSum(tempVelocities, &result));

    PetscScalar vsquared = result;

    result = (result - (this->numParticles - this->numConstraints) * this->temperature) / firstNoseHooverMass;
/*---------------------------------------------------------------------------------------------------------*/




    // PetscPrintf(PETSC_COMM_WORLD, "\nas1: %f, sum(V.^2): %f, N: %i, N_c: %i, T: %f, Q1: %f", result, vsquared, numParticles, numConstraints, temperature, firstNoseHooverMass);
/*---------------------------------------------------------------------------------------------------------*/

    /**
     * @todo This should be a Vec Copy?
    */
    PetscCall(VecSetValue(*accelerations, 0, result, INSERT_VALUES));

    PetscCall(VecAssemblyBegin(*accelerations));
    PetscCall(VecAssemblyEnd(*accelerations));

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

        /* Get the reduced values for each vector */
        const PetscScalar *evenNoseHooverMassArray;
        const PetscScalar *oddNoseHooverMassArray;
        const PetscScalar *noseHooverVelocitiesArray;

        PetscScalar *desiredMassesArray;
        PetscScalar *decrementedMassesArray;
        PetscScalar *decrementedVelocitiesArray;

        PetscCall(VecGetArrayRead(this->oddNoseHooverMass,  &oddNoseHooverMassArray));
        PetscCall(VecGetArrayRead(this->evenNoseHooverMass, &evenNoseHooverMassArray));
        PetscCall(VecGetArrayRead(noseHooverVelocities,     &noseHooverVelocitiesArray));

        PetscCall(VecGetArray(desiredMasses,            &desiredMassesArray));
        PetscCall(VecGetArray(decrementedMasses,        &decrementedMassesArray));
        PetscCall(VecGetArray(decrementedVelocities,    &decrementedVelocitiesArray));

        for (PetscInt i = 0; i < numReducedParticles; i++)
        {
            desiredMassesArray[i]           = oddNoseHooverMassArray[i+1];
            decrementedMassesArray[i]       = evenNoseHooverMassArray[i];
            decrementedVelocitiesArray[i]   = noseHooverVelocitiesArray[i];
        }

        /* Restore all values */
        PetscCall(VecRestoreArrayRead(this->oddNoseHooverMass,  &oddNoseHooverMassArray));
        PetscCall(VecRestoreArrayRead(this->evenNoseHooverMass, &evenNoseHooverMassArray));
        PetscCall(VecRestoreArrayRead(noseHooverVelocities,     &noseHooverVelocitiesArray));

        PetscCall(VecRestoreArray(desiredMasses,            &desiredMassesArray));
        PetscCall(VecRestoreArray(decrementedMasses,        &decrementedMassesArray));
        PetscCall(VecRestoreArray(decrementedVelocities,    &decrementedVelocitiesArray));
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
        /* Skip the first particle for odd results */
        PetscScalar* resultsArray;
        const PetscScalar* tempResultsArray;
        PetscInt numReducedParticles = this->NHChainOrder/2 - 1;

        PetscCall(VecGetArray(*result, &resultsArray));
        PetscCall(VecGetArrayRead(tempResults, &tempResultsArray));

        for (PetscInt i = 0; i < numReducedParticles; i++)
        {
            resultsArray[i+1] = tempResultsArray[i];
        }

        PetscCall(VecRestoreArray(*result, &resultsArray));
        PetscCall(VecRestoreArrayRead(tempResults, &tempResultsArray));
    }

    return errorStatus;
}


PetscErrorCode Hyperoptimization::calculatePositionIncrement(Vec previousPosition, Vec previousVelocity, PetscScalar timeStepIn, Vec *result)
{
    PetscErrorCode errorStatus = 0;

/*---------------------------------------------------------------------------------------------------------*/
        // PetscScalar firstPosition;
        // PetscScalar firstVelocity;
        // PetscScalar firstRight;
        // PetscScalar firstFinal;
        // PetscInt indexA = 0;
        // PetscCall(VecGetValues(previousPosition, 1, &indexA, &firstPosition));
        // PetscCall(VecGetValues(previousVelocity, 1, &indexA, &firstVelocity));
/*---------------------------------------------------------------------------------------------------------*/

    // Create a temporary vector for inner-results
    Vec tempVector;
    PetscCall(VecDuplicate(previousVelocity, &(tempVector)));
    PetscCall(VecCopy(previousVelocity, tempVector));

    /* dt * v[t-1] */
    PetscCall(VecScale(tempVector, timeStepIn));
/*---------------------------------------------------------------------------------------------------------*/
        // PetscCall(VecGetValues(tempVector, 1, &indexA, &firstRight));
/*---------------------------------------------------------------------------------------------------------*/

    /* x + (dt * v[t-1]) */
    PetscCall(VecAYPX(tempVector, 1, previousPosition));

    /* Copy return value */
    PetscCall(VecCopy(tempVector, *result));
/*---------------------------------------------------------------------------------------------------------*/
    // PetscCall(VecGetValues(tempVector, 1, &indexA, &firstFinal));
    // PetscPrintf(PETSC_COMM_WORLD, "\nPos: %f, Vel: %f, Right: %f, Final: %f, timeStepIn: %f", firstPosition, firstVelocity, meanPos);
/*---------------------------------------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------------------------------------*/
    // PetscScalar firstLeft;
    // PetscScalar firstRight;
    // PetscScalar firstFinal;
    // PetscInt indexA = 0;
    // PetscCall(VecGetValues(leftSide, 1, &indexA, &firstLeft));
    // PetscCall(VecGetValues(rightSide, 1, &indexA, &firstRight));
/*---------------------------------------------------------------------------------------------------------*/


    /* Add left and right and copy to result */
    PetscCall(VecAYPX(leftSide, 1, rightSide));
    PetscCall(VecCopy(leftSide, *result));

/*---------------------------------------------------------------------------------------------------------*/
    // PetscCall(VecGetValues(leftSide, 1, &indexA, &firstFinal));

    // PetscPrintf(PETSC_COMM_WORLD, "\nLS1: %f, RS1: %f, Tot: %ftimestep: %f\n", firstLeft, firstRight, firstFinal, timeStepIn);
/*---------------------------------------------------------------------------------------------------------*/

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

    return errorStatus;
}

PetscErrorCode Hyperoptimization::assembleNewPositions(PetscScalar firstNoseHooverVelocity, PetscScalar *lagrangianMultiplier)
{
    PetscErrorCode errorStatus = 0;

/*---------------------------------------------------------------------------------------------------------*/
    PetscReal minPos;
    PetscReal maxPos;
    PetscScalar meanPos;

    PetscCall(VecMax(newPosition, NULL, &maxPos));
    PetscCall(VecMin(newPosition, NULL, &minPos));
    PetscCall(VecMean(newPosition, &meanPos));

    // PetscPrintf(PETSC_COMM_WORLD, "\nmaxX: %f, minX: %f, meanX: %f", maxPos, minPos, meanPos);
/*---------------------------------------------------------------------------------------------------------*/

    /* Setup */
    Vec rightSide;
    PetscCall(VecDuplicate(this->newPosition, &(rightSide)));
    PetscCall(VecCopy(this->constraintSensitivities, rightSide));

    PetscInt numPositions;
    PetscCall(VecGetSize(this->newPosition, &numPositions));
    PetscInt indices[numPositions];
    for (PetscInt i = 0; i <= numPositions; i++)
    {
        indices[i] = i;
    }

    /* Lagrangian Multipliers require properly constrained design variables */
    PetscScalar *newPositionsArray;
    PetscCall(VecGetArray(this->newPosition, &newPositionsArray));
    for (PetscInt i = 0; i <= numPositions; i++)
    {
        if (newPositionsArray[i] > 1)
        {
            newPositionsArray[i] = 1;
        }
        if (newPositionsArray[i] < 0)
        {
            newPositionsArray[i] = 0;
        }
    }
    PetscCall(VecRestoreArray(this->newPosition, &newPositionsArray));

    // Fix all passive elements
    // opt->SetVariables(this->newPosition, opt->xPassive);


    /** @todo implement lagrangian multiplier */
    PetscScalar scaleFactor = - this->timestep * this->halfTimestep * std::exp(-this->halfTimestep * firstNoseHooverVelocity);
    PetscCall(VecScale(rightSide, scaleFactor));

/*---------------------------------------------------------------------------------------------------------*/
    // PetscReal minPos;
    // PetscReal maxPos;
    // PetscScalar meanPos;

    // PetscCall(VecMax(constraintSensitivities, NULL, &maxPos));
    // PetscCall(VecMin(constraintSensitivities, NULL, &minPos));
    // PetscCall(VecMean(constraintSensitivities, &meanPos));


    // PetscPrintf(PETSC_COMM_WORLD, "\ndt: %f, dgMax: %f, dgMin: %f, dgMean: %f, NHV1: %f\n", timestep, maxPos, minPos, meanPos, firstNoseHooverVelocity);

    // PetscCall(VecMax(newPosition, NULL, &maxPos));
    // PetscCall(VecMin(newPosition, NULL, &minPos));
    // PetscCall(VecMean(newPosition, &meanPos));

    // PetscPrintf(PETSC_COMM_WORLD, "\nmaxC: %f, minC: %f, meanC: %f", maxPos, minPos, meanPos);
/*---------------------------------------------------------------------------------------------------------*/

    PetscScalar lagrangianMult;
    this->lagMult.computeLagrangianMultiplier(this->newPosition, rightSide, this->numParticles, &lagrangianMult);
    PetscCall(VecAXPY(this->newPosition, lagrangianMult, rightSide));

    *lagrangianMultiplier = lagrangianMult;

    // Ensure new positions are all physical (restrict between 0 and 1)
    PetscCall(VecGetArray(this->newPosition, &newPositionsArray));

    for (PetscInt i = 0; i <= numPositions; i++)
    {
        if (newPositionsArray[i] > 1)
        {
            newPositionsArray[i] = 1;
        }
        if (newPositionsArray[i] < 0)
        {
            newPositionsArray[i] = 0;
        }
    }
    PetscCall(VecRestoreArray(this->newPosition, &newPositionsArray));

    // Fix all passive elements
    // opt->SetVariables(this->newPosition, opt->xPassive);

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateSensitvities(Vec positions)
{
    PetscErrorCode errorStatus = 0;
    // Positions used in constraint calculations
    PetscCall(VecCopy(positions, opt->x));
    errorStatus = this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);
    // CHKERRQ(errorStatus);

    // Compute sensitivities
    errorStatus = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), 
                                                                    &(opt->gx[0]),
                                                                    opt->dfdx,
                                                                    opt->dgdx[0],
                                                                    opt->xPhys,
                                                                    // opt->xPhys,
                                                                    opt->Emin,
                                                                    opt->Emax,
                                                                    opt->penal,
                                                                    opt->volfrac);//,
                                                                    // data);
    CHKERRQ(errorStatus);

    // Scale??

    // // Compute objective scale
    // if (firstIteration) {
    //     opt->fscale = 10.0 / opt->fx;
    //     firstIteration = false;
    // }

    // // Scale Sensitivities
    // VecScale(opt->dfdx, opt->fscale);

    // Calculate g and dgdx for the local volume constraint ??? 
    /** @todo ASK HAZHIR!! */
    // if (opt->localVolumeStatus) {
    //     ierr = local->Constraint(opt->xPhys, &(opt->gx[1]), opt->dgdx[1]);
    //     CHKERRQ(ierr);
    // }

    // Filter sensitivities (chainrule)
    errorStatus = filter->Gradients(opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, opt->beta,
                            opt->eta);
    CHKERRQ(errorStatus);

    PetscCall(VecCopy(opt->dfdx, this->sensitivities));
    PetscCall(VecScale(this->sensitivities, -1));
    PetscCall(VecCopy(opt->dgdx[0], this->constraintSensitivities));

    return errorStatus;
}

PetscErrorCode Hyperoptimization::calculateTemperature(Vec velocities, PetscScalar *temperature)
{
    PetscErrorCode errorStatus = 0;

    Vec tempVector;
    PetscCall(VecDuplicate(velocities, &tempVector));
    PetscCall(VecCopy(velocities, tempVector));

    PetscCall(VecPointwiseMult(tempVector, tempVector, tempVector));
    PetscCall(VecMean(tempVector, temperature));

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

    /* Setup all vectors */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(tempOddNoseHooverVelocity)));
    PetscCall(VecSetSizes(tempOddNoseHooverVelocity, PETSC_DETERMINE, numOddNHIndices));
    PetscCall(VecSetFromOptions(tempOddNoseHooverVelocity));

    PetscCall(VecDuplicate(this->newPosition, &(tempPosition)));
    // PetscCall(VecDuplicate(this->newPosition, &(prevPosition)));
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

    /* Initial Values */
    // PetscCall(VecCopy(this->positions.at(0), prevPosition));
    PetscCall(VecSet(prevEvenNoseHooverPosition,    0));
    PetscCall(VecSet(prevEvenNoseHooverVelocity,    0));
    PetscCall(VecSet(prevOddNoseHooverPosition,     0));
    PetscCall(VecSet(prevOddNoseHooverVelocity,     0));

    PetscPrintf(PETSC_COMM_WORLD, "...initialized!\n");

    for (PetscInt iteration = 0; iteration < this->numIterations; iteration++)
    {
        /* Calculate half-timestep positions */
/*---------------------------------------------------------------------------------------------------------*/
// PetscPrintf(PETSC_COMM_WORLD, "\n------------------------------------------------------------------\n");
/*---------------------------------------------------------------------------------------------------------*/
        calculatePositionIncrement(prevPosition, this->prevVelocity, this->halfTimestep, &tempPosition);
        calculatePositionIncrement(prevEvenNoseHooverPosition, prevEvenNoseHooverVelocity, this->halfTimestep, &tempEvenNoseHooverPosition);

        /* Calcualte previous Nose Hoover accelerations */
        calculateFirstNoseHooverAcceleration(this->prevVelocity, &tempOddNoseHooverAccelerations);
        calculateRemainingNoseHooverAccelerations(prevEvenNoseHooverVelocity, false, &tempOddNoseHooverAccelerations);

        /* Calculate half-timestep velocities */

/*---------------------------------------------------------------------------------------------------------*/
        // PetscReal minPos;
        // PetscReal maxPos;
        // PetscScalar meanPos;

        // PetscScalar firstNoseHooverA;
        // PetscInt indexA = 0;
        // PetscCall(VecGetValues(tempOddNoseHooverAccelerations, 1, &indexA, &firstNoseHooverA));

        // PetscCall(VecMax(tempOddNoseHooverAccelerations, NULL, &maxPos));
        // PetscCall(VecMin(tempOddNoseHooverAccelerations, NULL, &minPos));
        // PetscCall(VecMean(tempOddNoseHooverAccelerations, &meanPos));


        // PetscPrintf(PETSC_COMM_WORLD, "\nMaxPoNHA: %f, MinPoNHA: %f, MeanPoNHA: %f, firstNHA: %f", maxPos, minPos, meanPos, firstNoseHooverA);
/*---------------------------------------------------------------------------------------------------------*/

        calculateVelocityIncrement(prevOddNoseHooverVelocity, prevEvenNoseHooverVelocity, tempOddNoseHooverAccelerations, this->halfTimestep, &tempOddNoseHooverVelocity);


/*---------------------------------------------------------------------------------------------------------*/
        // PetscReal minTPos;
        // PetscReal maxTPos;
        // PetscScalar meanTPos;

        // PetscCall(VecMax(tempPosition, NULL, &maxTPos));
        // PetscCall(VecMin(tempPosition, NULL, &minTPos));
        // PetscCall(VecMean(tempPosition, &meanTPos));

        // PetscPrintf(PETSC_COMM_WORLD, "\nmaxTX: %f, minTX: %f, meanTX %f", maxTPos, minTPos, meanTPos);
/*---------------------------------------------------------------------------------------------------------*/

        /* Get objective and constraint sensitivities */
        calculateSensitvities(tempPosition);

        /* Get first nose hoover velocity */
        PetscScalar firstNoseHooverVelocity;
        const PetscScalar *tempOddNoseHooverVelocityArray;

        PetscCall(VecGetArrayRead(tempOddNoseHooverVelocity, &tempOddNoseHooverVelocityArray));
        firstNoseHooverVelocity = tempOddNoseHooverVelocityArray[0];
        PetscCall(VecRestoreArrayRead(tempOddNoseHooverVelocity, &tempOddNoseHooverVelocityArray));

        /* Calculate full increment values */
        calculateVelocityIncrement(this->prevVelocity, firstNoseHooverVelocity, this->sensitivities, this->timestep, &newVelocity);

/*---------------------------------------------------------------------------------------------------------*/
        // PetscReal minTPos;
        // PetscReal maxTPos;
        // PetscScalar meanTPos;

        // PetscCall(VecMax(sensitivities, NULL, &maxTPos));
        // PetscCall(VecMin(sensitivities, NULL, &minTPos));
        // PetscCall(VecMean(sensitivities, &meanTPos));

        // PetscReal minNV;
        // PetscReal maxNV;
        // PetscScalar meanNV;

        // PetscPrintf(PETSC_COMM_WORLD, "\nmaxSEN: %f, minSEN: %f, meanSEN: %f, NV: %f", maxTPos, minTPos, meanTPos, firstNoseHooverVelocity);
/*---------------------------------------------------------------------------------------------------------*/

        calculatePositionIncrement(prevOddNoseHooverPosition, tempOddNoseHooverVelocity, this->timestep, &newOddNoseHooverPosition);

        /* Nose Hoover acceleration at half-timestep */
        calculateRemainingNoseHooverAccelerations(tempOddNoseHooverVelocity, true, &tempEvenNoseHooverAccelerations);

        /* Add 0 to the end of the odd nose hoover velocities */
        Vec extendedOddNoseHooverVelocity;
        // const PetscScalar *tempOddNoseHooverVelocityArray;
        PetscScalar *extendedOddNoseHooverVelocityArray;

        PetscCall(VecDuplicate(tempOddNoseHooverVelocity, &extendedOddNoseHooverVelocity));
        PetscCall(VecGetArrayRead(tempOddNoseHooverVelocity, &tempOddNoseHooverVelocityArray));
        PetscCall(VecGetArray(extendedOddNoseHooverVelocity, &extendedOddNoseHooverVelocityArray));

        for (PetscInt i = 0; i < numOddNHIndices-1; i++)
        {
            extendedOddNoseHooverVelocityArray[i] = tempOddNoseHooverVelocityArray[i+1];
        }
        extendedOddNoseHooverVelocityArray[numOddNHIndices-1] = 0;

        // PetscCall(VecSetValue(extendedOddNoseHooverVelocity, numOddNHIndices-1, 0, INSERT_VALUES));
        // PetscCall(VecSetValues(extendedOddNoseHooverVelocity, numOddNHIndices-1, straightNHIndices, &(tempOddNoseHooverVelocityArray[1]), INSERT_VALUES));

        PetscCall(VecRestoreArrayRead(tempOddNoseHooverVelocity, &tempOddNoseHooverVelocityArray));
        PetscCall(VecRestoreArray(extendedOddNoseHooverVelocity, &extendedOddNoseHooverVelocityArray));
        // PetscCall(VecAssemblyBegin(extendedOddNoseHooverVelocity));
        // PetscCall(VecAssemblyEnd(extendedOddNoseHooverVelocity));

        /* Continue calculating full-timestep values */
        calculateVelocityIncrement(prevEvenNoseHooverVelocity, extendedOddNoseHooverVelocity, tempEvenNoseHooverAccelerations, this->timestep, &(newEvenNoseHooverVelocity));
        calculatePositionIncrement(tempPosition, newVelocity, this->halfTimestep, &(this->newPosition));

/*---------------------------------------------------------------------------------------------------------*/
        // PetscReal minTPos;
        // PetscReal maxTPos;
        // PetscScalar meanTPos;

        // PetscCall(VecMax(tempPosition, NULL, &maxTPos));
        // PetscCall(VecMin(tempPosition, NULL, &minTPos));
        // PetscCall(VecMean(tempPosition, &meanTPos));

        // PetscReal minNV;
        // PetscReal maxNV;
        // PetscScalar meanNV;

        // PetscCall(VecMax(newVelocity, NULL, &maxNV));
        // PetscCall(VecMin(newVelocity, NULL, &minNV));
        // PetscCall(VecMean(newVelocity, &meanNV));

        // PetscPrintf(PETSC_COMM_WORLD, "\nmaxTX: %f, minTX: %f, meanTX: %f, maxNV: %f, minNV: %f, meanNV: %f", maxTPos, minTPos, meanTPos, maxNV, minNV, meanNV);
/*---------------------------------------------------------------------------------------------------------*/

        calculatePositionIncrement(tempEvenNoseHooverPosition, newEvenNoseHooverPosition, this->halfTimestep, &(newEvenNoseHooverPosition));

        /* Calculate the new positions */
        PetscScalar lagMultiplier;
        assembleNewPositions(firstNoseHooverVelocity, &lagMultiplier);

        /* Calcualte the new velocities as v_new = (x_new -  x_old) / timestep */
        PetscCall(VecWAXPY(newVelocity, -1, prevPosition, newPosition));
        PetscCall(VecScale(newVelocity, 1/this->timestep));

        /* Calculate full-timestep Nose Hoover accelerations */
        calculateFirstNoseHooverAcceleration(newVelocity, &tempOddNoseHooverAccelerations);
        calculateRemainingNoseHooverAccelerations(newEvenNoseHooverPosition, false, &tempOddNoseHooverAccelerations);

        /* Calculate final odd Nose Hoover Velocity */
        calculateVelocityIncrement(tempOddNoseHooverVelocity, newEvenNoseHooverPosition, tempOddNoseHooverAccelerations, this->halfTimestep, &newOddNoseHooverVelocity);


        if (temperatureCheck)
        {
            PetscScalar temperature;
            calculateTemperature(prevVelocity, &temperature);

            PetscReal minPos;
            PetscReal maxPos;
            PetscScalar meanPos;

            PetscCall(VecMax(prevPosition, NULL, &maxPos));
            PetscCall(VecMin(prevPosition, NULL, &minPos));
            PetscCall(VecMean(prevPosition, &meanPos));

            PetscReal minVel;
            PetscReal maxVel;
            PetscScalar meanVel;

            PetscCall(VecMax(prevVelocity, NULL, &maxVel));
            PetscCall(VecMin(prevVelocity, NULL, &minVel));
            PetscCall(VecMean(prevVelocity, &meanVel));

            PetscPrintf(PETSC_COMM_WORLD, "iter: %i, Num Part: %i, Max Pos: %f, Min Pos: %f, Mean Pos: %f, Max Vel: %f, Min Vel: %f, Mean Vel: %f, Temperature: %f, LM: %f\n", iteration, numParticles, maxPos, minPos, meanPos, maxVel, minVel, meanVel, temperature, lagMultiplier);
        }

        /* SAVE VARIABLES! */
        // PetscCall(VecCopy(prevPosition, this->positions.at(iteration + 1)));

        PetscCall(VecCopy(newPosition, prevPosition));
        PetscCall(VecCopy(newVelocity, prevVelocity));

        PetscCall(VecCopy(newEvenNoseHooverPosition,    prevEvenNoseHooverPosition));
        PetscCall(VecCopy(newEvenNoseHooverVelocity,    prevEvenNoseHooverVelocity));
        PetscCall(VecCopy(newOddNoseHooverPosition,     prevOddNoseHooverPosition));
        PetscCall(VecCopy(newOddNoseHooverVelocity,     prevOddNoseHooverVelocity));
    }

    return errorStatus;
}
