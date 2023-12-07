#!/usr/bin/env python3

import HypOptLib

solver = HypOptLib.HypOptLib()

solver.setSavePath("test_timestepping.h5")

timestep = 0.0003

solver.setTargetTemperature(0.01)
solver.setTimestep(timestep)
solver.setMaximumIterations(2000)

solver.setRandomStartingValues(False)

# solver.enableVariableTimestep(1.04,          # timestepConstantAlpha
#                               0.98,         # timestepConstantBeta
#                               0.00000001)   # diffusionConstant

solver.loadInitialConditionsFromFile("../tests/randomInitial32x16x16_T0.01.h5")

solver.newRun(  [0,0],              # iterationSaveRange
                [32,16,16])         # gridDimensions

# solver.generateRandomInitialConditionsFile([32, 16, 16], "../tests/randomInitial32x16x16_T0.01.h5")

# solver.newRun(  False,              # randomStartingValues
#                 False,              # saveHamiltonian
#                 0.001,              # initialTemperature
#                 3.0,                # penalty
#                 0.08,               # minimumFilterRadius
#                 0.12,               # volumeFraction
#                 0.001,              # timestep
#                 10,                 # noseHooverChainOrder
#                 30,                 # maximumIterations
#                 [2,8],              # iterationSaveRange
#                 [32,16,16])         # gridDimensions

# solver.restartRun("testing.h5",
#                   10,
#                   [0,10],
#                   False)
