
#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include <petsc.h>
#include "TopOpt.h"
#include "LinearElasticity.h"
#include "Filter.h"
#include "HypOptException.h"
#include "Hyperoptimization.h"
#include "FileManager.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class HypOptLib
{
    public:
        HypOptLib(){}

        void setSavePath(std::string filePath)
        {
            this->savePath = filePath;
        }

        uint32_t newRun(bool                    randomStartingValues,
                        bool                    saveHamiltonian,
                        double                  initialTemperature,
                        double                  penalty,
                        double                  minimumFilterRadius,
                        double                  volumeFraction,
                        double                  timestep,
                        uint32_t                noseHooverChainOrder,
                        uint32_t                maximumIterations,
                        std::vector<uint32_t>  *iterationSaveRange,
                        std::vector<uint32_t>  *gridDimensions);

        uint32_t restartRun(std::string filePath,
                        uint32_t maximumIterations,
                        std::vector<uint32_t> *iterationSaveRange,
                        bool saveHamiltonian);

    private:
        PetscErrorCode runLoop(Hyperoptimization solver, PetscInt numItr, FileManager output);

        std::string savePath = "";
};


PYBIND11_MODULE(HypOptLib, m) 
{
    py::class_<HypOptLib>(m, "HypOptLib")
        .def(py::init<>())

        .def("newRun", &HypOptLib::newRun, "Initializes all parameters and starts a run")

        .def("restartRun",   &HypOptLib::restartRun, "Restartrs a Design Loop from a given file")

        .def("setSavePath", &HypOptLib::setSavePath, "Sets save path for all runs.");

    py::register_exception<HypOptException>(m, "HypOptError");
}
