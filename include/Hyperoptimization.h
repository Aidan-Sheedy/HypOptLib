/**
 * @file Hyperoptimization.h
 *
 * This file describes the layout of the Hyperoptimization class.
 *
 * @todo LICENSE
**/

#pragma once

#include "LinearElasticity.h"
#include "TopOpt.h"
#include "Filter.h"
// #include "topoptlib.h"

#include <petsc.h>
#include <vector>

class LagrangianMultiplier
{
    public:
        LagrangianMultiplier(){}

        LagrangianMultiplier(Filter *filter, TopOpt *opt)
        {
            this->filter    = filter;
            this->opt       = opt;
        }

        PetscScalar computeLagrangianMultiplier(Vec positions, Vec C, PetscInt numParticles, PetscScalar *returnValue)
        {
            PetscErrorCode errorStatus = 0;

            // Copies of input parameters
            Vec positionCopy;
            Vec CCopy;

            PetscCall(VecDuplicate(positions, &positionCopy));
            PetscCall(VecDuplicate(C, &CCopy));

            PetscCall(VecCopy(positions, positionCopy));
            PetscCall(VecCopy(C, CCopy));

            // Filter positions
            PetscCall(VecCopy(positionCopy, opt->x));
            // PetscCall(this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, opt->projectionFilter, opt->beta, opt->eta));
            PetscCall(this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta));
            PetscCall(VecCopy(opt->xPhys, positionCopy));

            // Filter positions
            PetscCall(VecCopy(CCopy, opt->x));
            // PetscCall(this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, opt->projectionFilter, opt->beta, opt->eta));
            PetscCall(this->filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta));
            PetscCall(VecCopy(opt->xPhys, CCopy));

            // Sum values
            PetscScalar positionSum;
            PetscCall(VecSum(positionCopy, &positionSum));
            PetscScalar CSum;
            PetscCall(VecSum(CCopy, &CSum));

            // PetscPrintf(PETSC_COMM_WORLD, "\nN: %i, volfrac: %f, posSum: %f, cSum: %f\n", numParticles, opt->volfrac, positionSum, CSum);

            *returnValue = (numParticles * opt->volfrac - positionSum )/ CSum;

            return errorStatus;
        }

    private:
        Filter *filter;

        TopOpt *opt;

};

/**
 * @class Hyperoptimization
 *
*/
class Hyperoptimization
{
    public:
        // Hyperoptimization(  LinearElasticity* physics,
        //                     TopOpt* opt,
        //                     Filter* filter,
        //                     DataObj data,
        //                     LagrangianMultiplier lagMult,
        //                     PetscScalar temperature,
        //                     Vec initialPositions,
        //                     Vec initialVelocities,
        //                     PetscScalar NHChainOrder,
        //                     PetscInt numIterations,
        //                     PetscFloat timestep);
        Hyperoptimization(){}

        PetscErrorCode init(LinearElasticity* physics,
                            TopOpt* opt,
                            Filter* filter,
                            // DataObj data,
                            LagrangianMultiplier lagMult,
                            PetscScalar temperature,
                            Vec initialPositions,
                            PetscScalar NHChainOrder,
                            PetscInt numIterations,
                            PetscFloat timestep);


        ~Hyperoptimization();

        PetscErrorCode runDesignLoop(); /* @note Maybe this should have some of the parameters? like initial values and such? Unsure. */


    private:
        /**
         * the acceleration of the first Nose Hoover particle. Solves the equation:
         * 
         * a_s1 = (sum(s_v1.^2) - (N - N_c) * T) / Q1
         * 
         * where:
         *  a_s1 =
         *  s_v1 =
         *  N    =
         *  N_c  =
         *  T    =
         *  Q1   =
         * 
         * @param allNoseHooverVelocities (Vec)
         * vector of Nose Hoover velocities with which to calculate the accelerations
         * 
         * @param accelerations (Vec*)
         * vector of accelerations in which to return the calculated acceleration (at the first index)
         *
         * @returns the calculated acceleration of the first Nose Hoover particle given the
         *          velocities provided (placed in the first index of accelerations).
         * 
        **/
        PetscErrorCode calculateFirstNoseHooverAcceleration(Vec allNoseHooverVelocities, Vec *accelerations);

        /**
         * Computes the acceleration of the all Nose Hoover particles beyond the first. Solves the
         * equation:
         * 
         * a_s{i} = (Q{i-1} * v_s{i-1}^2 - T)/Q{i}
         * 
         * where:
         *  a_s{i} = the i'th Nose Hoover acceleration
         *  v_s{i} = the i'th Nose Hoover velocity
         *  Q{i}   = the i'th Nose Hoover mass
         *  T      = the temperature
         * 
         * @param noseHooverVelocities (Vec)
         * v_s{i-1}, Vector of the Nose Hoover velocities. This is either the odd or even velocities
         * depending on the flag evenOutput.
         * 
         * @param evenOutput (bool)
         * Flag to indicate if even or odd acclerations are being calculated.
         * 
         * @param result (Vec*)
         * Array in which to add the return result.
         * 
         * @returns the accelerations calculated from the provided velocities. Returned in the pointer 
         *          to 'result'.
         *
        **/
        PetscErrorCode calculateRemainingNoseHooverAccelerations(Vec noseHooverVelocities, bool evenOutput, Vec *result);

        /**
         * Calculates particle positions for the provided timestep. Computes the equation:
         * 
         * x{t+dt} = x{t} + dt * v{t}
         * 
         * @param previousPosition (Vec)
         * Vector of positions at the previous timestep
         * 
         * @param previousVelocity (Vec)
         * Vector of velocities at the previous timestep
         * 
         * @param timeStep (PetscScalar)
         * The time at which to calculate the next position
         * 
         * @param result (Vec*)
         * Array in which to add the return result.
         * 
         * @returns the incremented position using the described equation. Returned to 
         *          the pointer to 'result'
        */
        PetscErrorCode calculatePositionIncrement(Vec previousPosition, Vec previousVelocity, PetscScalar timeStep, Vec *result);

        /**
         * Calculates the next velocity at the given time step. Completes the equation:
         * 
         * v_1{t+dt} = v_1{t} * exp(- dt * v_2) + dt * a * exp(- dt/2 * v_2)
         * 
         * @param velocityOne (Vec)
         * v_1{t}, the velocity to be incremented.
         * 
         * @param velocityTwo (Vec)
         * v_2, the second velocity. Typically even velocities if the incremented are odd,
         * and vice-versa.
         * 
         * @param acceleration (Vec)
         * a, the acceleration of the velocity to be incremented.
         * 
         * @param timeStep (PetscScalar)
         * dt, the timestep at which to calculate the new velocity
         * 
         * @param result (Vec*)
         * Array in which to add the return result.
         * 
         * @returns the incremented velocity using the described equation. Returned to 
         *          the pointer to 'result'
         * 
        **/
        PetscErrorCode calculateVelocityIncrement(Vec velocityOne, Vec velocityTwo, Vec acceleration, PetscScalar timeStep, Vec *result);

        PetscErrorCode calculateVelocityIncrement(Vec velocityOne, PetscScalar velocityTwo, Vec acceleration, PetscScalar timeStep, Vec *result);

        PetscErrorCode assembleNewPositions(PetscScalar firstNoseHooverVelocity, PetscScalar *lagrangianMultiplier);

        PetscErrorCode calculateTemperature(Vec velocities, PetscScalar *temperature);

        PetscErrorCode saveIteration(Vec positions, Vec velocities, Vec noseHooverPositions, Vec noseHooverVelocities, PetscScalar lagrangianMultiplier);

        PetscErrorCode calculateSensitvities(Vec positions);

    private:
        LinearElasticity* physics;

        TopOpt* opt;

        Filter* filter;

        // DataObj data;

        LagrangianMultiplier lagMult;

        PetscScalar temperature;

        PetscScalar NHChainOrder;

        PetscInt numParticles;

        PetscInt numIterations;

        PetscInt numConstraints;

        /** @todo FIGURE OUT HOW TO DO THIS */
        // std::vector<Vec> positions;

        /** @todo confirm if this is needed*/
        // std::vector<Vec> velocities;

        PetscFloat timestep;

        PetscFloat halfTimestep;

        /* Position of Nose Hoover particles */
        /** @todo confirm if this is needed*/
        // Vec evenNoseHooverPositions; 

        /* Velocities of Nose Hoover particles */
        /** @todo confirm if this is needed*/
        // Vec evenNoseHooverVelocities;

        /**
         * Mass of Nose Hoover particles.
         * 
         * @note    This is initialized automatically to 1 for now, 
         *          a better constructor should be made later.
        */
        Vec noseHooverMass;

        Vec evenNoseHooverMass;

        Vec oddNoseHooverMass;

        /** @todo confirm if this is needed*/
        // Vec passiveElements;

        Vec prevPosition;

        Vec prevVelocity;

        Vec newPosition;

        // Vec newVelocity;

        Vec sensitivities;

        Vec constraintSensitivities;

        /** @todo make this a pass-in variable/debugging parameter! */
        bool temperatureCheck = true;

        bool firstIteration = true;

        // Vec prevEvenNoseHooverPosition;

        // Vec prevEvenNoseHooverVelocity;

        // Vec newEvenNoseHooverPosition;

        // Vec newEvenNoseHooverVelocity;

        // Vec prevOddNoseHooverPosition;

        // Vec prevOddNoseHooverVelocity;

        // Vec newOddNoseHooverPosition;

        // Vec newOddNoseHooverVelocity;

};
