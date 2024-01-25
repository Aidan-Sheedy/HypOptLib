#!/usr/bin/env python3

import HypOptLib

solver = HypOptLib.HypOptLib()

solver.setSavePath("medium_run_5deg_restart_2.h5")

temperature = 5
timestep    = 0.00004
numItr      = 3600

solver.setTargetTemperature(temperature)
solver.setTimestep(timestep)
solver.setMaximumIterations(numItr)
solver.setVolumeFraction(0.12)
solver.setSaveFrequency(10)

# solver.setRandomStartingValues(False)
# solver.setMaxSimulationTime(2350)
# solver.setMaxSimulationTime(0.8)

# alpha = 1.4
# beta = 0.98

alpha = 1.115
beta = 0.9

# tempDifusionConst = 0.0000000000005
# volfracDiffusionConst = 0.0000001
volfracDiffusionConst   = 0.00000004

# solver.enableVariableTimestep(alpha,
#                               beta,
#                             #   tempDifusionConst)   # diffusionConstant
#                               volfracDiffusionConst)   # diffusionConstant -- volfrac version

# solver.setSaveHamiltonian(True)

# solver.loadInitialConditionsFromFile("../tests/randomInitial64x32x32_T5.h5")

# solver.newRun(  [0,10000],          # iterationSaveRange
#                 [64,32,32])         # gridDimensions

solver.restartRun("medium/t5/medium_run_5deg_restart_1.h5", [0,3600])

# solver.generateRandomInitialConditionsFile([64, 32, 32], "../tests/randomInitial64x32x32_T5.h5")
