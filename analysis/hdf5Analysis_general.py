import h5py
import numpy as np
import plotting
import matplotlib.pyplot as plt


def hdf5Analysis(filePaths, attributeKey, attributeName):

    dataGroup = "/Dataset"
    stateGroup = dataGroup + "/State"
    settingsGroup = "/Setting"

    outputFiles         = []
    settings            = []
    desiredTemperatures = []
    desiredAttributes   = []

    for i in range(0, len(filePaths)):
        outputFiles.append(h5py.File(filePaths[i]))
        settings.append(outputFiles[i][settingsGroup])

        desiredTemperatures.append(settings[i].attrs["T"])

        desiredAttributes.append(np.asarray(outputFiles[i][dataGroup + "/" + attributeKey]))

    plotting.plotAttributeSideBySide(desiredAttributes, desiredTemperatures, attributeName)

    plt.show() 

    return



filePath1 = '../run/testing_all_the_way.h5'

filePaths = [filePath1]

hdf5Analysis(filePaths, "Temperature", "Temperature")
