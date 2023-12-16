import h5py
import numpy as np
import matplotlib.pyplot as plt

def getHistogram(fileName, initialConditionsFile=False):
    file = h5py.File(fileName)
    if initialConditionsFile:
        finalPositions = file["Initial Position"]
        finalVelocities = file["Initial Velocity"]
    else:
        finalVelocities = file["Dataset/Final Velocity"]
        finalPositions = file["Dataset/Final Position"]

        computeTimes = np.asarray(file["Dataset/Iteration Compute Time"])
        timesteps = np.asarray(file["Dataset/Timestep"])

        totalSimTime = sum(timesteps)

        if 0 == totalSimTime:
            totalSimTime = file["Setting"].attrs['dt'] * file["Setting"].attrs['sampleN']

        print("--------------------------------------------------------------------------------------------")
        print("Stats for file:",fileName)
        print("Total compute time for simulation:", computeTimes.sum(), "Seconds, =", computeTimes.sum()/60, "Minutes = ", computeTimes.sum()/60/60, "hours.")
        print("Total simulation time:", totalSimTime)
        print("Average compute time/second of simTime:", computeTimes.sum()/totalSimTime, "\n")

    velocityCounts, velocityBins = np.histogram(finalVelocities, bins=100)
    positionCounts, positionBins = np.histogram(finalPositions, bins=100)

    return velocityCounts, velocityBins, positionCounts, positionBins


def initialConditionsHistograms(fileName):
    velocityCounts, velocityBins, positionCounts, positionBins = getHistogram(fileName, True)

    plt.stairs(velocityCounts, velocityBins)
    plt.title("Velocity Histogram")
    plt.show()

    plt.stairs(positionCounts, positionBins)
    plt.title("Position Histogram")
    plt.show()


def singleFileHistograms(fileName):
    velocityCounts, velocityBins, positionCounts, positionBins = getHistogram(fileName)

    plt.stairs(velocityCounts, velocityBins)
    plt.title("Velocity Histogram")
    plt.show()

    plt.stairs(positionCounts, positionBins)
    plt.title("Position Histogram")
    plt.show()

def compareHistograms(fileNames, title="", subtitles=[]):
    numFiles = len(fileNames)
    velocityCounts = []
    velocityBins = []
    positionCounts = []
    positionBins = []
    for i in range(0, numFiles):
        velCnt, velBin, posCnt, posBin = getHistogram(fileNames[i])
        velocityCounts.append(velCnt)
        velocityBins.append(velBin)
        positionCounts.append(posCnt)
        positionBins.append(posBin)
    
    figure, subplots = plt.subplots(numFiles, sharex=True)
    figure.suptitle(title + " - Veloctiy\nNot at exact same time, not directly comparable.")
    figure.supylabel("Count")
    figure.supxlabel("Velocity")

    for i in range(0, numFiles):
        subplots[i].stairs(velocityCounts[i], velocityBins[i])
        subplots[i].set_ylim([0, np.sort(velocityCounts[i])[-2]])
        if (len(subtitles)==len(fileNames)):
            subplots[i].set_title(subtitles[i])
    plt.show()

    figure, subplots = plt.subplots(numFiles, sharex=True)
    figure.suptitle(title + " - Position\nNot at exact same time, not directly comparable.")
    figure.supylabel("Count")
    figure.supxlabel("Position")
    for i in range(0, numFiles):
        subplots[i].stairs(positionCounts[i], positionBins[i])
        subplots[i].set_ylim([0, np.sort(positionCounts[i])[-2]])
        if (len(subtitles)==len(fileNames)):
            subplots[i].set_title(subtitles[i])
    plt.show()


# ##################################################################################################
# # T = 0
# compareHistograms(["../run/timestep_tests/t0/test_timestepping_temp_vartime_t0.h5",
#                    "../run/timestep_tests/t0/test_timestepping_volfrac_vartime_t0.h5",
#                    "../run/timestep_tests/t0/test_timestepping_straightime_t0.h5"],
#                    "Temp=0",
#                    ["temp dt", "volfrac dt", "const dt"])

##################################################################################################
# T = 0.0001
# compareHistograms(["../run/timestep_tests/t0.0001/test_timestepping_temp_vartime_t0.0001.h5",
#                    "../run/timestep_tests/t0.0001/test_timestepping_volfrac_temptime_t0.0001.h5",
#                    "../run/timestep_tests/t0.0001/test_timestepping_straightime_t0.0001.h5"],
#                    "Temp=0.0001",
#                    ["temp dt", "volfrac dt", "const dt"])
    
##################################################################################################
# T = 1
compareHistograms(["../run/timestep_tests/t1/test_timestepping_volfrac_vartime_t1.h5",
                #    "../run/timestep_tests/t0/test_timestepping_volfrac_vartime_t0.h5",
                   "../run/test_timestepping_straightime_t1 (1).h5"],
                   "Temp=1",
                   ["volfrac dt", "const dt"])
