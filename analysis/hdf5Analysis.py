

import h5py
import numpy as np
import plotting
import matplotlib.pyplot as plt

# filePath = '../outputs/hypopt_output_first_run_broke.h5'
filePath = '../bin/hypopt_output.h5'
dataGroup = "/Dataset"
stateGroup = dataGroup + "/State"
settingsGroup = "/Setting"

outputFile = h5py.File(filePath)

settings = outputFile[settingsGroup]
# data = outputFile[dataGroup]
timestep = settings.attrs["dt"]
desiredTemperature = settings.attrs["T"]

temperature     = outputFile[dataGroup + "/Temperature"]
lagrangians     = outputFile[dataGroup + "/Lambda"]
should_be_vf    = outputFile[dataGroup + "/should_be_vf"]
should_be_vf_pt = outputFile[dataGroup + "/should_be_vf post truncate"]

plotting.plotAttribute(timestep, np.asarray(temperature), desiredTemperature, "temperature", ylimUpper=20)
plotting.plotAttribute(timestep, np.asarray(lagrangians), desiredTemperature, "Lagrangian Multipliers", ylimUpper=0.1e8, ylimLower=-1e5)
plotting.plotAttribute(timestep, np.asarray(should_be_vf), desiredTemperature, "should_be_vf", ylimUpper=0.120005, ylimLower=0.1199999)
plotting.plotAttribute(timestep, np.asarray(should_be_vf_pt), desiredTemperature, "should_be_vf post truncate", ylimUpper=0.120001, ylimLower=0.11995)

plt.show()

def findMaxMin(iteration, stateGroup):
    iterationName = stateGroup + "/iteration" + str(iteration)
    state = outputFile[iterationName]

    state = np.asarray(state)

    maxState = state.max()
    minState = state.min()

    print("Iteration ", iteration, " max: ", maxState, ", min: ", minState)

    return

# findMaxMin(800, stateGroup)
# for i in range(0, 3):
#     findMaxMin(i, stateGroup)