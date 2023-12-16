
// #include "HypOptLib.h"
// #include <vector>

// int main(int argc, char* argv[])
// {
//     HypOptLib solver;

//     solver.setSavePath("test_overflow.h5");

//     double timestep = 0.000005;

//     solver.setTargetTemperature(1);
//     solver.setTimestep(timestep);
//     solver.setMaximumIterations(10000);

//     solver.setRandomStartingValues(false);
//     solver.setMaxSimulationTime(1);

//     // # alpha = 1.4
//     // # beta = 0.98

//     double alpha = 1.115;
//     double beta = 0.9;

//     double tempDifusionConst = 0.0000000000005;
//     double volfracDiffusionConst = 0.000000001;

//     solver.loadInitialConditionsFromFile("../tests/randomInitial32x16x16_T1.h5");

//     std::vector<uint32_t> saveRange = {0,0};
//     std::vector<uint32_t> dimension = {32,16,16};

//     solver.newRun(&saveRange, &dimension);

// }
