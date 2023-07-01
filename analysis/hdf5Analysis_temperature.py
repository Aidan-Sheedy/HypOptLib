import h5py
import numpy as np
import plotting
import matplotlib.pyplot as plt

filePath1 = '../outputs/hypopt_output_first_run_broke.h5'
filePath2 = '../outputs/hypopt_output_smallscale_system.h5'
dataGroup = "/Dataset"
stateGroup = dataGroup + "/State"
settingsGroup = "/Setting"

outputFile1 = h5py.File(filePath1)
outputFile2 = h5py.File(filePath2)

settings1 = outputFile1[settingsGroup]
settings2 = outputFile2[settingsGroup]

timestep1 = settings1.attrs["dt"]
desiredTemperature1 = settings1.attrs["T"]
desiredTemperature2 = settings2.attrs["T"]

temperature1 = np.asarray(outputFile1[dataGroup + "/Temperature"])
temperature2 = np.asarray(outputFile2[dataGroup + "/Temperature"])

temps = [temperature1, temperature2]
desiredTemps = [desiredTemperature1, desiredTemperature2]

plotting.plotAttributeSideBySide(timestep1, temps, desiredTemps, "Temperature")

plt.show()
