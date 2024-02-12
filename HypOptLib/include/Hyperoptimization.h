/***************************************************************************//**
 * @file Hyperoptimization.h
 *
 * This file describes the layout of the Hyperoptimization class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#pragma once

#include <petsc.h>
#include <vector>
#include <limits>
#include "SensitivitiesWrapper.h"
#include "FilterWrapper.h"
#include "LagrangeMultiplier.h"
#include "HypOptParameters.h"
#include "FileManager.h"



/**
 * Hyperoptimization implementation.
 *
 * Main functionality is in initialization and design loop functions. This class assumes
 * position vectors are set up in accordance to the settings used by whatever Petsc solver is
 * being used.
 *
 * The class is abstracted from Petsc as much as possible, other than the use of Petsc type variables.
 * Wrappers are provided for all operations that are problem specific, i.e. for filtering, Lagrangian
 * multiplier calculation, and sensitivity calculation. This is so that while the implementation provided
 * is intended for topology optimization, these wrappers can be modified for any optimization problem that
 * uses Petsc.
 */
class Hyperoptimization
{
    public:
        /**
         * Empty constructor. Initialization functions are used to comply with Petsc error handling.
         */
        Hyperoptimization(){}

        /**
         * Destructor.
         *
         * @todo set this up properly.
         */
        ~Hyperoptimization()
        {
            VecDestroy(&evenNoseHooverMass);
            VecDestroy(&oddNoseHooverMass);
            VecDestroy(&newPosition);
            VecDestroy(&sensitivities);
            VecDestroy(&constraintSensitivities);

            VecDestroy(&prevState.evenNoseHooverPosition);
            VecDestroy(&prevState.evenNoseHooverVelocity);
            VecDestroy(&prevState.oddNoseHooverPosition);
            VecDestroy(&prevState.oddNoseHooverVelocity);
            VecDestroy(&prevState.position);
            VecDestroy(&prevState.velocity);
        }

        /**
         * Initialization funciton.
         *
         * Sets up the design loop. This version initializes the Nose Hoover chain and is intended
         * for new runs with no restart variables. Chains the initialization by calling the overloaded
         * initialization function.
         *
         * @param sensitivitiesWrapper wrapper called to calculate sensitivities.
         * @param filter wrapper called to calculate all filtering.
         * @param lagMult Lagrangian Multiplier wrapper.
         * @param temperature target temperature for the system.
         * @param initialPositions position vectors initialized to desired distribution.
         * @param initialVelocities velocity vectors initialized to desired distribution.
         * @param NHChainOrder number of Nose Hoover particles.
         * @param numIterations number of iterations to simulate.
         * @param timestep timestep between iterations.
         * @param fileManager contains helper functions for saving states and iterations.
         * @param iterationSaveRange range between which to save position data.
         * @param saveHamiltonian optional parameter to save Hamiltonian and Compliance data.
         * @param volumeFraction target volume fraction.
         * @param saveFrequency frequency at which to save position state within save range.
         * @param maxSimTime if set, maxmimum amount of simulation time to iterate for.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode init(SensitivitiesWrapper& sensitivitiesWrapper,
                            FilterWrapper& filter,
                            LagrangeMultiplier& lagMult,
                            PetscScalar temperature,
                            Vec initialPositions,
                            Vec initialVelocities,
                            PetscInt NHChainOrder,
                            PetscInt numIterations,
                            PetscScalar timestep,
                            FileManager* fileManager,
                            std::vector<uint32_t> iterationSaveRange,
                            bool saveHamiltonian,
                            PetscScalar volumeFraction,
                            uint32_t saveFrequency,
                            double maxSimTime = std::numeric_limits<double>::max());

        /**
         * Overloaded initialization funciton.
         *
         * Sets up the design loop. This version takes restarted Nose Hoover chain variables.
         *
         * @param sensitivitiesWrapper wrapper called to calculate sensitivities.
         * @param lagMult Lagrangian Multiplier wrapper.
         * @param filter wrapper called to calculate all filtering.
         * @param fileManager contains helper functions for saving states and iterations.
         * @param NHChainOrder number of Nose Hoover particles.
         * @param timestep timestep between iterations.
         * @param temperature target temperature for the system.
         * @param numIterations number of iterations to simulate.
         * @param iterationSaveRange range between which to save position data.
         * @param initialPositions position vectors initialized to desired distribution.
         * @param initialVelocities velocity vectors initialized to desired distribution.
         * @param initialEvenNoseHooverPosition initial positions of Nose Hoover particles with even indices
         * @param initialEvenNoseHooverVelocity initial velocities of Nose Hoover particles with even indices
         * @param initialOddNoseHooverPosition initial positions of Nose Hoover particles with odd indices
         * @param initialOddNoseHooverVelocity initial velocities of Nose Hoover particles with odd indices
         * @param saveHamiltonian optional parameter to save Hamiltonian and Compliance data.
         * @param volumeFraction target volume fraction.
         * @param saveFrequency frequency at which to save position state within save range.
         * @param maxSimTime if set, maxmimum amount of simulation time to iterate for.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode init(SensitivitiesWrapper& sensitivitiesWrapper,
                            LagrangeMultiplier& lagMult,
                            FilterWrapper& filter,
                            FileManager* fileManager,
                            PetscInt NHChainOrder,
                            PetscScalar timestep,
                            PetscScalar temperature,
                            PetscInt numIterations,
                            std::vector<uint32_t> iterationSaveRange,
                            Vec initialPositions,
                            Vec initialVelocities,
                            Vec initialEvenNoseHooverPosition,
                            Vec initialEvenNoseHooverVelocity,
                            Vec initialOddNoseHooverPosition,
                            Vec initialOddNoseHooverVelocity,
                            bool saveHamiltonian,
                            PetscScalar volumeFraction,
                            uint32_t saveFrequency,
                            double maxSimTime = std::numeric_limits<double>::max());

        /**
         * Sets up variable timestepping.
         *
         * @todo fill out info on how the algorithm works when Hazhir is ready.
         *
         * @param timestepConstantAlpha
         * @param timestepConstantBeta
         * @param timestepConstantK
         * @param diffusionConstant
         */
        void enableVariableTimestep(PetscScalar timestepConstantAlpha,
                                    PetscScalar timestepConstantBeta,
                                    PetscScalar diffusionConstant);

        /**
         * Main iteration loop which implements the hyperoptimization design algorithm.
         *
         * @todo basic info on the loop here?
         */
        PetscErrorCode runDesignLoop();

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @todo throw exception if save Hamiltonian is false?
         *
         * @returns vector of Hamiltonians for each iteration.
         */
        std::vector<PetscScalar> getHamiltonians() {return hamiltonians;}

        /**
         * Accessor for final vector of Compliance after running design loop.
         *
         * @todo throw exception if save Hamiltonian is false?
         *
         * @returns vector of Compliance for each iteration.
         */
        std::vector<PetscScalar> getCompliance() {return compliance;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns vector of temperatures for each iteration.
         */
        std::vector<PetscScalar> getTemperatures() {return temperatures;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns vector of Lagrange Multipliers for each iteration.
         */
        std::vector<PetscScalar> getLagrangeMultipliers() {return LagrangeMultipliers;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns vector of iteration times for each iteration.
         */
        std::vector<PetscScalar> getIterationTimes() {return iterationTimes;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns vector of iteration times for each iteration.
         */
        std::vector<PetscInt> getSolverIterationsHamiltonian() {return hamiltonianSolverIterations;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns vector of iteration times for each iteration.
         */
        std::vector<PetscInt> getSolverIterationsSensitivity() {return sensitivitySolverIterations;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns vector of positions for the final iteration.
         */
        HypOptParameters getFinalState() {return prevState;}

        /**
         * Accessor for final vector of Hamiltonians after running design loop.
         *
         * @returns save hamiltonian option.
         */
        bool getSaveHamiltonian() {return saveHamiltonian;}

        /**
         * Accessor for timestep value at each iteration.
         *
         * @returns a vector of the timestep used at each iteration.
         */
        std::vector<PetscScalar> getTimesteps() {return timesteps;}

        /**
         * Accessor for energy error at each timestep.
         *
         * @todo update description when variable timestep is finished.
         *
         * @returns a vector of the energy error at each iteration.
         */
        std::vector<PetscScalar> getEnergyErrors() {return energyErrors;}

        /**
         * Accessor for volume fraction at each iteration.
         *
         * @returns a vector of the systems volume fraction calculated at each iteration
         */
        std::vector<PetscScalar> getVolFracs() {return volfracs;}

    private:
        /**
         * the acceleration of the first Nose Hoover particle. Solves the equation:
         *
         * @f[
         * a_{s_1} = \frac{\sum_i(v_i^2) - (N - N_c)T}{Q_1}
         * @f]
         *
         * where:
         *  - \f$ a_{s1} \f$ is the first Nose Hoover particle's acceleration
         *  - \f$ v_i    \f$ is the i'th design particle's velocity
         *  - \f$ N      \f$ is number of particles
         *  - \f$ N_c    \f$ is number of constraints
         *  - \f$ T      \f$ is target temperature
         *  - \f$ Q_1    \f$ is mass of the first Nose Hoover particle
         *
         * @param allVelocities vector of particle velocities with which to calculate the accelerations
         * @param accelerations [out] vector of Nose Hoover accelerations (only first index populated)
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateFirstNoseHooverAcceleration(Vec allVelocities, Vec *accelerations);

        /**
         * Computes the acceleration of the all Nose Hoover particles beyond the first. Solves the
         * equation:
         *
         * @f[
         * a_{s_i} = \frac{Q_{i-1} v_{s_{i-1}}^2 - T}{Q_i}
         * @f]
         *
         * where:
         *  - \f$ a_{s_i} \f$ is the i'th Nose Hoover acceleration
         *  - \f$ v_{s_i} \f$ is the i'th Nose Hoover velocity
         *  - \f$ Q{i}    \f$ is the i'th Nose Hoover mass
         *  - \f$ T       \f$ is the target temperature
         *
         * @param noseHooverVelocities \f$ v_{s_{i-1}} \f$, Vector of the Nose Hoover velocities. This is either the odd or even velocities depending on the flag evenOutput.
         * @param evenOutput flag to indicate if even or odd acclerations are being calculated.
         * @param result [out] resulting accelerations. They are populated in the correct locations in the array. I.e. for odd output, the positions
         *               populated are \f$ i = 3, 5, 7, ... \f$.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateRemainingNoseHooverAccelerations(Vec noseHooverVelocities, bool evenOutput, Vec *result);

        /**
         * Calculates particle positions for the provided timestep. Computes the equation:
         *
         * @f[
         * x(t + \Delta t) = x(t_a) + \Delta t  v(t_b),
         * @f]
         *
         * where:
         *  - \f$ x(t + \Delta t) \f$ is the incremented position,
         *  - \f$ x(t_a) \f$ is the position at a previous time \f$ t_a \f$,
         *  - \f$ \Delta t \f$ is the timestep,
         *  - \f$ v(t_b) \f$ is the velocity at a previous time \f$ t_b \f$.
         *
         * where the positions and velocities can be for either design particles or Nose Hoover particles.
         * Note that \f$ t_a \f$ and \f$ t_b \f$ are typically either the same or off by a half timestep.
         *
         * @param previousPosition \f$ x(t_a) \f$
         * @param previousVelocity \f$ v(t_b) \f$
         * @param timeStep \f$ \Delta t \f$
         * @param result [out] array in which to add the return resulting position vector. Must match shape of input vectors.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculatePositionIncrement(Vec previousPosition, Vec previousVelocity, PetscScalar timeStep, Vec *result);

        /**
         * Calculates the next velocity at the given time step. Completes the equation:
         *
         * @f[
         * v_1(t + \Delta t) = v_1(t) e^{-\Delta t  v_2(t_a)} + \Delta t a(t_a) e^{-\frac{\Delta t}{2} v_2(t_a)},
         * @f]
         *
         * where:
         *  - \f$ v_1(t + \Delta t) \f$ is the incremented primary velocity,
         *  - \f$ v_1(t) \f$ is the previous primary velocity,
         *  - \f$ \Delta t \f$ is the timestep
         *  - \f$ v_2(t_a) \f$ is a secondary velocity term (odd Nose Hoover if the primary is even, and vice versa)
         *                     at a time \f$ t_a \f$
         *  - \f$ a(t_a) \f$ is the primary velcities' first derivative at the time \f$ t_a \f$.
         *
         * @param velocityOne \f$ v_1(t) \f$, vector of the velocity to be incremented.
         * @param velocityTwo \f$ v_2(t_a) \f$, vector of the second velocities. Typically even velocities if the incremented are odd, and vice-versa.
         * @param acceleration \f$ a(t_a) \f$, vector of the acceleration of the velocity to be incremented.
         * @param timeStep \f$ \Delta t \f$, the timestep at which to calculate the new velocity.
         * @param result [out] array in which to add the return result. Must match shape of input vectors.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateVelocityIncrement(Vec velocityOne, Vec velocityTwo, Vec acceleration, PetscScalar timeStep, Vec *result);

        /**
         * Overloaded function to calculates the next velocity at the given time step.
         *
         * The only difference from the base function is that \f$ v_2(t_a) \f$ is allowed to be a scalar value. This is necessary as
         * one step requires using the first Nose Hoover velocity for this term for every element in the resulting vector.
         *
         * @see Hyperoptimization::calculateVelocityIncrement(Vec velocityOne, Vec velocityTwo, Vec acceleration, PetscScalar timeStep, Vec *result)
         *
         * @param velocityOne \f$ v_1(t) \f$, vector of the velocity to be incremented.
         * @param velocityTwo \f$ v_2(t_a) \f$, scalar second velocity.
         * @param acceleration \f$ a(t_a) \f$, vector of the acceleration of the velocity to be incremented.
         * @param timeStep \f$ \Delta t \f$, the timestep at which to calculate the new velocity.
         * @param result [out] array in which to add the return result. Must match shape of input vectors.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateVelocityIncrement(Vec velocityOne, PetscScalar velocityTwo, Vec acceleration, PetscScalar timeStep, Vec *result);

        /**
         * Calcualtes the final incremented position vector. This step is somewhat complicated, but can be broken down as follows:
         *
         * 1. Truncate the positions at \f$ t + \frac{\Delta t}{2} \f$ to be within the range \f$ [0,1] \f$.
         * 2. Scale the constraint sensitivities by a calculated amount.
         * 3. Calculate Lagrangian multiplier.
         * 4. Scale the step 2 result by the Lagrangian multiplier
         * 5. Pointwise add the new position to the result of step 4
         * 6. Truncate the result of step 5 to get a meaningful value to the new positions.
         *
         * The final value for the new positions is saved in the parameter Hyperoptimization#newPosition.
         *
         * @note The full derivation of this step is more complicated, and the hyperoptimization paper should be consulted for more details.
         *
         * @param firstNoseHooverVelocity the first Nose Hoover velocity at time \f$ t + \frac{\Delta t}{2} \f$
         * @param LagrangeMultiplier [out] the Lagrangian multiplier is calculated as a byproduct of this function, and is returned here. Can be NULL if not desired.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode assembleNewPositions(PetscScalar firstNoseHooverVelocity, PetscScalar *LagrangeMultiplier);

        /**
         * Calcualtes the system temperature given the vector of design particle velocities.
         *
         * The temperature is calculated as:
         *
         * @f[
         * T = \frac{\sum_{i=0}^{N}{v_i^2}}{N}
         * @f]
         *
         * where:
         *  - \f$ T \f$ is the system temperature.
         *  - \f$ N \f$ is the number of design particles.
         *  - \f$ v_i \f$ is the i'th design particle velocity.
         *
         * @param velocities the design particle velocities from which to calculate temperature.
         * @param temperature [out] resulting temperature to populate.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateTemperature(Vec velocities, PetscScalar *temperature);

        /**
         * Calculates the Hamiltonian and objective function of the system.
         *
         * The Hamiltonian is calculated as:
         *
         * @f[
         * H = O + \frac{v^2}{2},
         * @f]
         *
         * where:
         *  - \f$ H \f$ is the Hamiltonian,
         *  - \f$ O \f$ is the objective function, which is abstracted and calculated by the problem-specific SensitivitiesWrapper,
         *  - \f$ v \f$ is the vector of disgn particle velocities.
         *
         * @param velocities vector of disgn particle velocities.
         * @param positions vector of disgn particle positions (used for calculating the objective function).
         * @param hamiltonian [out] resulting Hamiltonian to populate.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateHamiltonian(Vec velocities, Vec positions, PetscScalar *hamiltonian);

        /**
         * Calculates both constraint and objective sensitivities.
         *
         * The provided design particle positions are first filtered, before calculating the sensitivities using the
         * abstracted-sensitivity wrapper. Finally, the sensitivities are filtered. The results are stored in
         * Hyperoptimization::sensitivities and Hyperoptimization::constraintSensitivities respectively.
         *
         * @param positions vector of disgn particle positions.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode calculateSensitvities(Vec positions);

        /**
         * Truncates the design particle positions to be real-valued, within the range [0,1].
         *
         * All values greater than 1 are truncated to exactly 1, and all values less than 0 are truncated to
         * exactly 0. While this does add some small error to the resulting volume fraction and velocity distributions,
         * it is sufficiently small as to not affect the design iteration process. It is interesting to note that low
         * temperature systems will tend to have many particles close to zero, meaning that more particles are truncated
         * from negative values to 0 than from high values to 1. This gives a high tendancy to slightly err high on the
         * objective function.
         */
        PetscErrorCode truncatePositions(Vec *positions);

        /**
         * Calculates the next timestep for the variable timestep algorithm.
         *
         * Uses the following equation:
         *
         * @f[
         * \Delta t_i = \alpha \beta^k \Delta t_{i=i}
         * @f]
         *
         * where:
         *  - \f$ \Delta t_i \f$ is the timestep at step \f$ i \f$
         */
        void calculateNextTimeStep();

    private:
        FileManager* fileManager;
        SensitivitiesWrapper sensitivitiesWrapper; /** @todo find a better name than sensitivitiesWrapper */
        FilterWrapper filter;
        LagrangeMultiplier lagMult;
        PetscScalar temperature;
        PetscInt NHChainOrder;
        PetscInt numParticles;
        PetscInt numIterations;
        PetscInt numConstraints;
        PetscScalar timestep;
        PetscScalar halfTimestep;
        Vec evenNoseHooverMass; /** Mass of Nose Hoover particles.  @note This is initialized automatically to 1 for now, a better constructor should be made later. */
        Vec oddNoseHooverMass;
        HypOptParameters prevState;
        Vec newPosition;
        Vec sensitivities;
        Vec constraintSensitivities;
        std::vector<PetscScalar> LagrangeMultipliers;
        std::vector<PetscScalar> hamiltonians;
        std::vector<PetscScalar> compliance;
        std::vector<PetscScalar> temperatures;
        std::vector<PetscScalar> iterationTimes;
        std::vector<PetscInt> sensitivitySolverIterations;
        std::vector<PetscInt> hamiltonianSolverIterations;
        verbosity printInfo = INFO; /** @todo make this a pass-in variable/debugging parameter! */
        bool doneSolving = false;
        bool saveHamiltonian;

        bool saveRangeUseSimTime = false;
        uint32_t saveFrequency;
        std::vector<uint32_t> iterationSaveRange;

        bool variableTimestep = false;
        PetscScalar previousTimestep = 0;
        PetscScalar timestepConstantAlpha = 1.1;
        PetscScalar timestepConstantBeta = 0.99;
        PetscScalar timestepConstantK = 1;
        PetscScalar diffusionConstant = 0.00000001;
        PetscScalar previousTemperature = 0;
        double      maxSimTime;

        std::vector<PetscScalar> timesteps;
        std::vector<PetscScalar> energyErrors;
        std::vector<PetscScalar> volfracs;
        PetscScalar volumeFraction;
};
