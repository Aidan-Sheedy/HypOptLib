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
#include "FilterWrapper.h"
#include "LagrangeMultiplier.h"

#include <petsc.h>
#include <vector>

/**
 * @class Hyperoptimization
 *
*/
class Hyperoptimization
{
    public:
        Hyperoptimization(){}

        PetscErrorCode init(LinearElasticity* physics,
                            TopOpt* opt,
                            FilterWrapper* filter,
                            // DataObj data,
                            LagrangeMultiplier lagMult,
                            PetscScalar temperature,
                            Vec initialPositions,
                            Vec initialVelocities,
                            PetscScalar NHChainOrder,
                            PetscInt numIterations,
                            PetscScalar timestep);

        PetscErrorCode init(LinearElasticity* physics,
                            TopOpt* opt,
                            FilterWrapper* filter,
                            // DataObj data,
                            LagrangeMultiplier lagMult,
                            PetscScalar temperature,
                            Vec initialPositions,
                            Vec initialVelocities,
                            PetscScalar NHChainOrder,
                            PetscInt numIterations,
                            PetscScalar timestep,
                            PetscInt numIterationsToSave,
                            bool saveHamiltonian);


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

        PetscErrorCode assembleNewPositions(PetscScalar firstNoseHooverVelocity, PetscScalar *LagrangeMultiplier);

        PetscErrorCode calculateTemperature(Vec velocities, PetscScalar *temperature);

        PetscErrorCode calculateHamiltonian(Vec velocities, Vec positions, PetscScalar *hamiltonian);

        PetscErrorCode initializeHDF5();

        PetscErrorCode saveIteration(PetscInt iteration, Vec positions);

        PetscErrorCode saveFinalValues();

        PetscErrorCode calculateSensitvities(Vec positions);

        PetscErrorCode truncatePositions(Vec *positions);

    private:
        LinearElasticity* physics;

        TopOpt* opt;

        FilterWrapper* filter;

        LagrangeMultiplier lagMult;

        std::string saveFilePath;

        const std::string stateGroup = "/Dataset/State";

        const std::string dataGroup = "/Dataset";

        PetscScalar temperature;

        PetscScalar NHChainOrder;

        PetscInt numParticles;

        PetscInt numIterations;

        PetscInt numConstraints;

        PetscScalar timestep;

        PetscScalar halfTimestep;

        /**
         * Mass of Nose Hoover particles.
         * 
         * @note    This is initialized automatically to 1 for now, 
         *          a better constructor should be made later.
        */
        Vec evenNoseHooverMass;

        Vec oddNoseHooverMass;

        /** @todo confirm if this is needed*/
        // Vec passiveElements;

        Vec prevPosition;

        Vec prevVelocity;

        Vec newPosition;

        Vec sensitivities;

        Vec constraintSensitivities;

        std::vector<PetscScalar> LagrangeMultipliers;

        std::vector<PetscScalar> hamiltonians;

        std::vector<PetscScalar> compliance;

        std::vector<PetscScalar> genericData2;

        std::vector<PetscScalar> genericData;

        std::vector<PetscScalar> temperatures;

        std::vector<PetscScalar> iterationTimes;

        PetscInt numIterationsToSave;

        /** @todo make this a pass-in variable/debugging parameter! */
        bool temperatureCheck = true;

        bool saveData = true;

        bool saveHamiltonian;
};
