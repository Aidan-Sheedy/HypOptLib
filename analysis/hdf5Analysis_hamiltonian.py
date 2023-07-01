import h5py
import numpy as np
import plotting
import matplotlib.pyplot as plt

filePath1 = '../outputs/hypopt_output_FILENAME_HERE.h5'
filePath2 = '../outputs/hypopt_output_FILENAME_HERE.h5'
filePath3 = '../outputs/hypopt_output_FILENAME_HERE.h5'
dataGroup = "/Dataset"
stateGroup = dataGroup + "/State"
settingsGroup = "/Setting"

outputFile1 = h5py.File(filePath1)
outputFile2 = h5py.File(filePath2)
outputFile3 = h5py.File(filePath3)

settings1 = outputFile1[settingsGroup]
settings2 = outputFile2[settingsGroup]
settings3 = outputFile3[settingsGroup]

timestep1 = settings1.attrs["dt"]
desiredTemperature1 = settings1.attrs["T"]
desiredTemperature2 = settings2.attrs["T"]
desiredTemperature3 = settings3.attrs["T"]

hamiltonian1 = np.asarray(outputFile1[dataGroup + "/Hamiltonian"])
hamiltonian2 = np.asarray(outputFile2[dataGroup + "/Hamiltonian"])
hamiltonian3 = np.asarray(outputFile3[dataGroup + "/Hamiltonian"])

hamiltonians = [hamiltonian1, hamiltonian2, hamiltonian3]
desiredTemps = [desiredTemperature1, desiredTemperature2, desiredTemperature3]

plotting.plotAttributeSideBySide(timestep1, hamiltonians, desiredTemps, "Hamiltonian")

plt.show()
