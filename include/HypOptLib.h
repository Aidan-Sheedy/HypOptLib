
#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include <petsc.h>
#include "TopOpt.h"
#include "LinearElasticity.h"
#include "Filter.h"
#include "HypOptException.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class HypOptLib
{
    public:
        HypOptLib(){}
        void init(  bool                    randomStartingValues,
                    bool                    saveHamiltonian,
                    double                  initialTemperature,
                    double                  penalty,
                    double                  minimumFilterRadius,
                    double                  volumeFraction,
                    double                  timestep,
                    uint32_t                noseHooverChainOrder,
                    uint32_t                maximumIterations,
                    std::vector<double>    *iterationSaveRange,
                    std::vector<uint32_t>  *gridDimensions);



        // void init(std::string restartFilePath);

        void startLoop();


    private:
        PetscErrorCode initNoRestart(bool randomStartingValues);

        std::string help = "3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n";

        TopOpt* opt;

        LinearElasticity* physics;

        Filter* filter;

        bool saveHamiltonian;

        double initialTemperature;

        double timestep;

        double volumeFraction;

        uint32_t noseHooverChainOrder;

        uint32_t maximumIterations;

        std::vector<double> iterationSaveRange;

        std::vector<uint32_t> gridDimensions;

        Vec initialPositions;

        Vec initialVelocities;

};


PYBIND11_MODULE(HypOptLib, m) 
{
    py::class_<HypOptLib>(m, "HypOptLib")
        .def(py::init<>())

        .def("init", &HypOptLib::init, "Initializes all parameters")

        .def("startLoop",   &HypOptLib::startLoop, "Begins Design Loop");

    py::register_exception<HypOptException>(m, "HypOptError");
}
