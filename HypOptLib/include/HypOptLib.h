/***************************************************************************//**
 * @file HypOptLib.h
 *
 * Contains the main HypOptLib class, along with the pybind11 wrapper.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
 ******************************************************************************/

#pragma once

#include <stdint.h>
#include <vector>
#include <string>
#include <limits>

#include <petsc.h>
#include "LinearElasticity.h"
#include "Filter.h"
#include "HypOptException.h"
#include "Hyperoptimization.h"
#include "FileManager.h"
#include "HypOptParameters.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * Main HypOptLib class.
 *
 * This is the equivalent to a "main" function, and
 * provides access to the Python wrapper to the essential functions to run
 * Hyperoptimization.
 */
class HypOptLib
{
    public:
        /**
         * An empty constructor is used to make Pybind11 translation easier.
         */
        HypOptLib(){}

        /**
         * Starts a fresh run with the provided parameters.
         *
         * @param iterationSaveRange    Range of iterations to save. Does not need to include the final iteration to support restarting, this is saved regardless.
         *
         * @throws HypOptException
         *
         * @return 0 on success, or error. Only used to align with Petsc exception handling.
         *
         * @todo Some of these settings might lend themselves better to optional parameters, such as "setSavePath".
         */
        uint32_t newRun(std::vector<uint32_t>  *iterationSaveRange);

        /**
         * Restarts a simualtion from the provided file path. All options are parsed from the metadata in the restart file.
         *
         * @param filePath Location of the restart file.
         * @param iterationSaveRange Range of iterations to save. Does not need to include the final iteration to support restarting, this is saved regardless.
         *
         * @throws HypOptException
         *
         * @return 0 on success, or error. Only used to align with Petsc exception handling.
         *
         * @todo Some of these settings might lend themselves better to optional parameters, such as "setSavePath".
         */
        uint32_t restartRun(std::string restartPath,
                            std::vector<uint32_t> *iterationSaveRange);

        /**
         * Generates a file with randomized initial velocities and positions. This file can then be passed
         * as an optional parameter to ensure the same initial conditions are used accross multiple runs.
         *
         * @param gridDimensions Dimensions of the cells in the grid. The final mesh will have be one higher in each dimension. The dimenions should each be divisible by 2 three times.
         * @param filePath Name and path to save the output to.
         *
         * @throws HypOptException
         */
        void generateRandomInitialConditionsFile(std::vector<uint32_t> *gridDimensions, std::string filePath);

        /**
         * Optional setting. Sets the save path and file name of the simulation result.
         *
         * @note Defaults to 'hypopt_output.h5'
         */
        void setSavePath(std::string filePath)
        {
            this->savePath = filePath;
        }

        /**
         * Recommended setting. Sets the target system temperature for the simulation on new runs. Must be greater than 0.
         *
         * @note defaults to 0.
         */
        void setTargetTemperature(double targetTemperature)
        {
            if (0 > targetTemperature)
            {
                throw HypOptException("Invalid target temperature, must be greater than 0.");
            }
            this->targetTemperature = targetTemperature;
        }

        /**
         * Recommended setting. Sets the timestep between each iteration.
         *
         * @note defaults to 0.001.
         */
        void setTimestep(double timestep)
        {
            if (0 > timestep)
            {
                throw HypOptException("Invalid timestep, must be greater than 0.");
            }
            this->timestep = timestep;
        }

        /**
         * Recommended setting. Sets the number of Nose Hoover particles in the chain.
         *
         * @todo confirm if this has to be even, if yes set a check here
         *
         * @note defaults to 10.
         */
        void setNoseHooverChainOrder(uint32_t noseHooverChainOrder)
        {
            this->noseHooverChainOrder = noseHooverChainOrder;
        }

        /**
         * Recommended setting. Sets the number of Nose Hoover particles in the chain.
         *
         * @todo topopt settings does NOT currently get this (silence it?)
         *
         * @note defaults to 100.
         */
        void setMaximumIterations(uint32_t maximumIterations)
        {
            this->maximumIterations = maximumIterations;
        }

        /**
         * Optional setting. Sets optimization penalty power.
         *
         * @note defaults to 3.
         */
        void setPenalty(double penalty)
        {
            this->penalty = penalty;
        }

        /**
         * Optional setting. Sets the minimum filter radius. Must be greater than 0.
         *
         * @note defaults to 0.08.
         */
        void setMinimumFilterRadius(double minimumFilterRadius)
        {
            if (0 > minimumFilterRadius)
            {
                throw HypOptException("Invalid minimum filter radius, must be greater than 0.");
            }
            this->minimumFilterRadius = minimumFilterRadius;
        }

        /**
         * Optional setting. Sets the volume fraction.
         *
         * @todo topopt settings does NOT currently get this value.
         *
         * @note defaults to 0.12.
         */
        void setVolumeFraction(double volumeFraction)
        {
            if (0 > volumeFraction)
            {
                throw HypOptException("Invalid volume fraction, must be greater than 0.");
            }
            this->volumeFraction = volumeFraction;
        }

        /**
         * Debugging setting. If false, all positions be initialized to the volume fraction.
         *
         * @note Defaults to true.
         */
        void setRandomStartingValues(bool randomStartingValues)
        {
            this->randomStartingValues = randomStartingValues;
        }

        /**
         * Debugging parameter. Enabling will double run time, but save the compliance and Hamiltonian.
         *
         * @note Defaults to false.
         */
        void setSaveHamiltonian(bool saveHamiltonian)
        {
            this->saveHamiltonian = saveHamiltonian;
        }

        /**
         * Sets the maximum amount of time the simulation can run for.
         *
         * This is calculated as timestep*iteration. Useful for variable timestep runs when the number of iterations might not
         * be predictable.
         *
         * @note If this is set, the lower of simulation time and setMaximumIterations will be used.
         *
         * @param simTime Maxmimum allowed simulation time.
         */
        void setMaxSimulationTime(double simTime)
        {
            this->maxSimTime = simTime;
        }

        /**
         * Sets the frequency that iterations are saved. Still only saved within the iteration save range.
         * 
         * For example, if set to 3 then the system state is only saved every three iterations.
         *
         * @note Defaults to 1, meaning every iteration is saved.
         *
         * @param frequency Frequency at which to save system states.
         */
        void setSaveFrequency(double frequency)
        {
            this->saveFrequency = frequency;
        }

        /**
         * Enables variable timestepping with the provided parameters.
         *
         * @param timestepConstantAlpha Growth parameter. Should be close to but greater than 1.
         * @param timestepConstantBeta Decay parameter. Should be close to but less than 1.
         * @param diffusionConstant Conditional error parameter.
         */
        void enableVariableTimestep(double timestepConstantAlpha,
                                    double timestepConstantBeta,
                                    double diffusionConstant)
        {
            this->timestepConstantAlpha = timestepConstantAlpha;
            this->timestepConstantBeta = timestepConstantBeta;
            this->diffusionConstant = diffusionConstant;
            this->variableTimestep = true;
        }

        /**
         * Loads the starting positions and velocities from the provided file.
         *
         * @param filePath HDF5 file containing a valid initial conditions given the temperature and grid dimensions provided.
         */
        void loadInitialConditionsFromFile(std::string filePath)
        {
            /** @todo check if this is a real path somewhere down the line */
            this->initialConditionsFile = filePath;
            this->initialConditionsFromFile = true;
            this->randomStartingValues = false;
        }

        /**
         * Sets the provided dimensions for any future simulations.
         *
         * @param gridDimensions x, y, and z coordinate dimensions for the grid. The mesh will be 1 larger in each direction.
         * @param domain domain coordinates for the grid. This defines the actual unit lengths of each dimension. For example, a
         * mesh with grid dimensions [32, 32, 32] and domain [(0, 1), (0, 1), (-1, 1)] means that each cell in the grid will be a
         * rectangle, with the z axis twice the width of other two dimensions.
         */
        void setGridProperties(std::vector<uint32_t> *gridDimensions, DomainCoordinates domain);

        /**
         * Sets the provided list of boundary conditions.
         *
         * Each boundary condition is defined by its:
         *  - Type: FIXED_POINT or LOAD.
         *  - <x,y,z>Range: range in each dimension over which the condition applies. This can either be a range, ie 
         *      x coordinates from [0, 4] or a single value, ie [4, 4]. Only points with x coordinate 4 will apply this BC.
         *  - degreesOfFreedom: a set of either 0, 1, or 2, where 0 is the x axis, 1 is the y axis, and 2 is the z axis.
         *  - value: for LOAD types, the scalar value for the load.
         *
         * @param boundaryConditions list of boundary conditions to apply. If an empty or invalid list is supplied, this will throw an error.
         */
        void setBoundaryConditions(std::vector<BoundaryCondition> *boundaryConditions);

        /**
         * Compatibility flag. Call this function if the restart file provided was generated before the custom mesh and boundary
         * conditions were implemented.
         *
         * @warning if this flag is set, the grid properties and boundary conditions supplied *must* match the restart file.
         */
        void restartDoesntSupportCustomMesh()
        {
            this->restartUseNewMeshes = true;
        }

    private:
        /**
         * Utility function, runs the design loop set up by either newRun or restartRun.
         *
         * @param solver used to run the design loop.
         * @param numItr number of iterations to run.
         * @param output HDF5 file manager.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode runLoop(Hyperoptimization solver, PetscInt numItr, FileManager output);

        /**
         * Randomizes the initial positions and velocities.
         *
         * @param initialPosition [out] array of initial positions to populate.
         * @param initialVelocity [out] array of initial velocities to populate.
         *
         * @todo explain how the algorithms work.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        PetscErrorCode randomizeStartingVectors(Vec initialPosition, Vec initialVelocity);

        std::string savePath                = "hypopt_output.h5";
        double      targetTemperature       = 0;
        double      penalty                 = 3;
        double      minimumFilterRadius     = 0.08;
        double      volumeFraction          = 0.12;
        double      timestep                = 0.001;
        uint32_t    noseHooverChainOrder    = 10;
        uint32_t    maximumIterations       = 100;
        uint32_t    saveFrequency           = 1;
        bool        randomStartingValues    = true;
        bool        saveHamiltonian         = false;
        bool        restartUseNewMeshes     = false;

        std::vector<uint32_t> gridDimensions;
        DomainCoordinates domain = {0};
        std::vector<BoundaryCondition> boundaryConditions;

        bool   variableTimestep = false;
        double timestepConstantAlpha = 1.1;
        double timestepConstantBeta = 0.99;
        double diffusionConstant = 0.00000001;
        double maxSimTime = std::numeric_limits<double>::max();

        bool        initialConditionsFromFile = false;
        std::string initialConditionsFile = "";
};

/**
 * Python wrapper. Defines functions which are directly accessible to the python library.
 *
 * The wrapper exposes the following class and functions:
 *
 * **class** *HypOptLib* class responsible for initializing parameters and starting simulations.
 *
 * **function** *newRun* Initializes all parameters and starts a run
 *
 * **function** *restartRun* Restarts a Design Loop from a given file
 *
 * **function** *setSavePath* Sets save path for all runs.
 *
 * **function** *setTargetTemperature* Optional setting. Sets the save path and file name of the simulation result.
 *
 * **function** *setPenalty* Debugging setting. If false, all positions be initialized to the volume fraction.
 *
 * **function** *setMinimumFilterRadius* Optional setting. Sets optimization penalty power.
 *
 * **function** *setVolumeFraction* Optional setting. Sets the minimum filter radius. Must be greater than 0.
 *
 * **function** *setTimestep* Recommended setting. Sets the timestep between each iteration.
 *
 * **function** *setNoseHooverChainOrder* Recommended setting. Sets the number of Nose Hoover particles in the chain.
 *
 * **function** *setMaximumIterations* Recommended setting. Sets the number of Nose Hoover particles in the chain.
 *
 * **function** *setRandomStartingValues* Debugging setting. If false, all positions be initialized to the volume fraction.
 *
 * **function** *setSaveHamiltonian* Debugging parameter. Enabling will double run time, but save the compliance and Hamiltonian.
 */
PYBIND11_MODULE(HypOptLib, m)
{
    py::class_<HypOptLib>(m, "HypOptLib")
        .def(py::init<>())

        .def("newRun",                  &HypOptLib::newRun,                 "Initializes all parameters and starts a run")
        .def("restartRun",              &HypOptLib::restartRun,             "Restarts a Design Loop from a given file")
        .def("setSavePath",             &HypOptLib::setSavePath,            "Sets save path for all runs.")
        .def("setTargetTemperature",    &HypOptLib::setTargetTemperature,   "Optional setting. Sets the save path and file name of the simulation result.")
        .def("setPenalty",              &HypOptLib::setPenalty,             "Debugging setting. If false, all positions be initialized to the volume fraction.")
        .def("setMinimumFilterRadius",  &HypOptLib::setMinimumFilterRadius, "Optional setting. Sets optimization penalty power.")
        .def("setVolumeFraction",       &HypOptLib::setVolumeFraction,      "Optional setting. Sets the minimum filter radius. Must be greater than 0.")
        .def("setTimestep",             &HypOptLib::setTimestep,            "Recommended setting. Sets the timestep between each iteration.")
        .def("setNoseHooverChainOrder", &HypOptLib::setNoseHooverChainOrder,"Recommended setting. Sets the number of Nose Hoover particles in the chain.")
        .def("setMaximumIterations",    &HypOptLib::setMaximumIterations,   "Recommended setting. Sets the number of Nose Hoover particles in the chain.")
        .def("setRandomStartingValues", &HypOptLib::setRandomStartingValues,"Debugging setting. If false, all positions be initialized to the volume fraction.")
        .def("setSaveHamiltonian",      &HypOptLib::setSaveHamiltonian,     "Debugging parameter. Enabling will double run time, but save the compliance and Hamiltonian.")
        .def("enableVariableTimestep",  &HypOptLib::enableVariableTimestep, "Enables variable timestepping with the provided parameters.")
        .def("generateRandomInitialConditionsFile",  &HypOptLib::generateRandomInitialConditionsFile, "Generates an HDF5 file with randomized initial position and velocity vectors.")
        .def("loadInitialConditionsFromFile",  &HypOptLib::loadInitialConditionsFromFile, "Optional setting to load initial conditions from a file.")
        .def("setSaveFrequency",        &HypOptLib::setSaveFrequency,       "Optional setting to change the frequency at which system states are saved within the iteration save range.")
        .def("setMaxSimulationTime",    &HypOptLib::setMaxSimulationTime,   "Optional setting to finish simulation at a given time.")
        .def("setGridProperties",       &HypOptLib::setGridProperties,      "Required settings for grid dimensions and domain.")
        .def("setBoundaryConditions",   &HypOptLib::setBoundaryConditions,  "Required setting for problem boundary conditions.")
        .def("restartDoesntSupportCustomMesh",   &HypOptLib::restartDoesntSupportCustomMesh,  "Optional setting for cross compatibility.");

    py::enum_<BoundaryConditionType>(m, "BoundaryConditionType")
        .value("FIXED_POINT",   FIXED_POINT)
        .value("LOAD",          LOAD)
        .export_values();

    py::class_<BoundaryCondition>(m, "BoundaryCondition")
        .def(py::init<>())
        .def_readwrite("type",      &BoundaryCondition::type)
        .def_readwrite("xRange",    &BoundaryCondition::xRange)
        .def_readwrite("yRange",    &BoundaryCondition::yRange)
        .def_readwrite("zRange",    &BoundaryCondition::zRange)
        .def_readwrite("degreesOfFreedom", &BoundaryCondition::degreesOfFreedom)
        .def_readwrite("value",     &BoundaryCondition::value);

    py::class_<DomainCoordinates>(m, "DomainCoordinates")
        .def(py::init<>())
        .def_readwrite("xMinimum", &DomainCoordinates::xMinimum)
        .def_readwrite("xMaximum", &DomainCoordinates::xMaximum)
        .def_readwrite("yMinimum", &DomainCoordinates::yMinimum)
        .def_readwrite("yMaximum", &DomainCoordinates::yMaximum)
        .def_readwrite("zMinimum", &DomainCoordinates::zMinimum)
        .def_readwrite("zMaximum", &DomainCoordinates::zMaximum);

    py::register_exception<HypOptException>(m, "HypOptError");
}
