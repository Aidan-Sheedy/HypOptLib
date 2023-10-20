#!/usr/bin/env python3

import HypOptLib

solver = HypOptLib.HypOptLib()

solver.setSavePath("testing_all_the_way.h5")

solver.newRun(  False,              # randomStartingValues
                False,              # saveHamiltonian
                0.001,              # initialTemperature
                3.0,                # penalty
                0.08,               # minimumFilterRadius
                0.12,               # volumeFraction
                0.001,              # timestep
                10,                 # noseHooverChainOrder
                30,                 # maximumIterations
                [2,8],              # iterationSaveRange
                [32,16,16])         # gridDimensions

# solver.restartRun("testing.h5",
#                   10,
#                   [0,10],
#                   False)
