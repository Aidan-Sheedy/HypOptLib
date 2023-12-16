import h5py
import numpy as np
import matplotlib.pyplot as plt



def compareAttribute(fileNames, attributeName, title="", subtitles=[]):
    numFiles = len(fileNames)

    timesteps = []
    simTimes = []
    attributes = []

    for fileName in fileNames:
        file = h5py.File(fileName)
        timestep = np.asarray(file["Dataset" + "/Timestep"])
        attribute = np.asarray(file["Dataset" + "/" + attributeName])

        if (0 == sum(timestep)):
            dt = file["Setting"].attrs["dt"]
            numItr = len(attribute)
            simTime = dt * numItr
            timestep = np.linspace(0, numItr*dt, numItr)
        else:
            simTime = sum(timestep)
            time = 0
            for i in range(0, len(attribute)):
                time += timestep[i]
                timestep[i] += time

        timesteps.append(timestep)
        simTimes.append(simTime)
        attributes.append(attribute)

    # Cut short at shortest simtime
    for i in range(0, numFiles):
        time = 0
        j = -1
        while(time < min(simTimes)):
            j += 1
            time = timesteps[i][j]
        # Ignore the first timestep as it is usually very large
        timesteps[i] = timesteps[i][1:j]
        attributes[i] = attributes[i][1:j]

    figure, subplots = plt.subplots(numFiles, sharex=True)
    figure.suptitle(title)
    figure.supylabel(attributeName)
    figure.supxlabel("Time")

    for i in range(0, numFiles):
        subplots[i].plot(timesteps[i], attributes[i], marker='.')
        if (len(subtitles)==len(fileNames)):
            subplots[i].set_title(subtitles[i])

    plt.show()

##############################################################################################
# t = 0
# compareAttribute(["../run/timestep_tests/t0/test_timestepping_temp_vartime_t0.h5",
#                    "../run/timestep_tests/t0/test_timestepping_volfrac_vartime_t0.h5",
#                    "../run/timestep_tests/t0/test_timestepping_straightime_t0.h5"],
#                    "Temperature",
#                    "Temp=0",
#                    ["temp dt", "volfrac dt", "const dt"])

##############################################################################################
# t = 0.0001
# variableTimestep = "../run/timestep_tests/t0.0001/test_timestepping_temp_vartime_t0.0001.h5"
# straight = "../run/timestep_tests/t0.0001/test_timestepping_straightime_t0.0001.h5"
# varTime2 = "../run/timestep_tests/t0.0001/test_timestepping_volfrac_temptime_t0.0001.h5"

# compareAttribute([variableTimestep, varTime2, straight], "Temperature", "Temp=0.0001", ["temp dt", "volfrac dt", "const dt"])


##############################################################################################
# t = 0.01
# compareAttribute(["../run/timestep_tests/t0.01/old/test_timestepping_tempvar_vartime_t0.01.h5",
#                   "../run/timestep_tests/t0.01/old/test_timestepping_volfrac_vartime_t0.01.h5",
#                   "../run/timestep_tests/t0.01/old/test_timestepping_straighttime_t0.01.h5"],
#                    "Temperature",
#                    "Temp=0.01",
#                    ["temp dt", "volfrac dt", "const dt"])

# compareAttribute(["../run/timestep_tests/t0.01/test_timestepping_volfrac_vartime_t0.01.h5",
#                   "../run/timestep_tests/t0.01/old/test_timestepping_volfrac_vartime_t0.01.h5"],
#                    "Temperature",
#                    "Temp=0.01",
#                    ["max_volfracerr=1e-7", "max_volfracerr=1e-10"])

# compareAttribute(["../run/timestep_tests/t0.01/test_timestepping_volfrac_vartime_t0.01.h5",
#                   "../run/test_timestepping_straightime_t0.01.h5"],
#                    "Temperature",
#                    "Temp=0.01",
#                    ["vartime volfrac", "straight"])

##############################################################################################
# t = 1
compareAttribute(["../run/timestep_tests/t1/test_timestepping_volfrac_vartime_t1.h5",
                #    "../run/timestep_tests/t0/test_timestepping_volfrac_vartime_t0.h5",
                   "../run/test_timestepping_straightime_t1 (1).h5"],
                   "Temperature",
                   "Temp=0",
                   ["temp dt", "volfrac dt", "const dt"])

