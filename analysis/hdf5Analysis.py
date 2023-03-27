import h5py
import numpy as np
import plotting
import matplotlib.pyplot as plt


def hdf5Analysis(filePaths, attributeKey, attributeName, distinguishingFeature="Temperature", title="", single=False, saveFig=False):

    dataGroup = "/Dataset"
    stateGroup = dataGroup + "/State"
    settingsGroup = "/Setting"

    outputFiles         = []
    settings            = []
    distinguishers      = []
    desiredAttributes   = []

    for i in range(0, len(filePaths)):
        outputFiles.append(h5py.File(filePaths[i]))
        settings.append(outputFiles[i][settingsGroup])

        if (distinguishingFeature == "Temperature"):
            distinguishers.append(settings[i].attrs["T"])
        elif(distinguishingFeature == "Particles"):
            volume = settings[i].attrs["nelx"] * settings[i].attrs["nely"] * settings[i].attrs["nelz"]
            distinguishers.append(volume)
        else:
            raise Exception("Invalid Distinguishing Feature: " + distinguishingFeature + 
                            "\nValid features are: `Temperature` `Particles`")

        desiredAttributes.append(np.asarray(outputFiles[i][dataGroup + "/" + attributeKey]))

    print("Plotting", title)

    if (single):
        plotting.plotAttribute(1, desiredAttributes[0], distinguishers[0], attributeKey)
    else:
        plotting.plotAttributeSideBySide(desiredAttributes, distinguishers, distinguishingFeature, attributeName, title, savePlots=saveFig)

    return


# filePath2 = '../outputs/hypopt_output_first_run_broke.h5'
# filePath1 = '../outputs/hypopt_output_smallscale_system.h5'

# filePaths = [filePath1, filePath2]

# hdf5Analysis(filePaths, "Temperature", "Temperature", "Particles", "System Temperature at 0 Degrees", saveFig=True)

filePath = ['../bin/hypopt_output.h5']#, '../bin/hypopt_output_NO_LAGRANGIAN_FILTERING.h5']
# filePath = ['../bin/hypopt_output_NO_LAGRANGIAN_FILTERING.h5']

hdf5Analysis(filePath, "Temperature",       "Temperature"           , single=True)
hdf5Analysis(filePath, "Lambda",            "Lagrangian Multiplier" , single=True)
hdf5Analysis(filePath, "Hamiltonian",       "Hamiltonian"           , single=True)
hdf5Analysis(filePath, "Compliance",        "Compliance"            , single=True)
hdf5Analysis(filePath, "Volume Fraction",   "Volume Fraction"       , single=True)
hdf5Analysis(filePath, "Max Position",      "Max Position"          , single=True)

plt.show()