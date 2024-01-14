#!/usr/bin/env python3

import HypOptLib

solver = HypOptLib.HypOptLib()

solver.setSavePath("test_timestepping_straightime_t1.h5")

timestep = 0.000005

solver.setTargetTemperature(1)
solver.setTimestep(timestep)
solver.setMaximumIterations(10)

solver.setRandomStartingValues(False)
solver.setMaxSimulationTime(10)

# alpha = 1.4
# beta = 0.98

alpha = 1.115
beta = 0.9

tempDifusionConst = 0.0000000000005
volfracDiffusionConst = 0.000000001

# solver.enableVariableTimestep(alpha,
#                               beta,
#                             #   tempDifusionConst)   # diffusionConstant
#                               volfracDiffusionConst)   # diffusionConstant -- volfrac version

# solver.loadInitialConditionsFromFile("../tests/randomInitial32x16x16_T1.h5")

solver.newRun(  [0,0],              # iterationSaveRange
                [32,16,16])         # gridDimensions

# solver.generateRandomInitialConditionsFile([32, 16, 16], "../tests/randomInitial32x16x16_T0.0001.h5")


def test():

    return