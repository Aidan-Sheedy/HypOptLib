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

#include <petsc.h>
#include "LinearElasticity.h"
#include "Filter.h"
#include "HypOptException.h"
#include "Hyperoptimization.h"
#include "FileManager.h"

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
         * @param gridDimensions        Dimension of cells in the grid. The final mesh will have be one higher in each dimension.
         *
         * @return 0 on success, or error. Only used to align with Petsc exception handling.
         *
         * @todo Some of these settings might lend themselves better to optional parameters, such as "setSavePath".
         */
        uint32_t newRun(std::vector<uint32_t>  *iterationSaveRange,
                        std::vector<uint32_t>  *gridDimensions);

        /**
         * Restarts a simualtion from the provided file path. All options are parsed from the metadata in the restart file.
         *
         * @param filePath Location of the restart file.
         * @param iterationSaveRange Range of iterations to save. Does not need to include the final iteration to support restarting, this is saved regardless.
         *
         * @return 0 on success, or error. Only used to align with Petsc exception handling.
         *
         * @todo Some of these settings might lend themselves better to optional parameters, such as "setSavePath".
         */
        uint32_t restartRun(std::string restartPath,
                            std::vector<uint32_t> *iterationSaveRange);

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

        std::string savePath                = "hypopt_output.h5";
        double      targetTemperature       = 0;
        double      penalty                 = 3;
        double      minimumFilterRadius     = 0.08;
        double      volumeFraction          = 0.12;
        double      timestep                = 0.001;
        uint32_t    noseHooverChainOrder    = 10;
        uint32_t    maximumIterations       = 100;
        bool        randomStartingValues    = true;
        bool        saveHamiltonian         = false;
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
        .def("setSaveHamiltonian",      &HypOptLib::setSaveHamiltonian,     "Debugging parameter. Enabling will double run time, but save the compliance and Hamiltonian.");

    py::register_exception<HypOptException>(m, "HypOptError");
}
