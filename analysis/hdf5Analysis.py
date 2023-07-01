import h5py
import numpy as np
import plotting
import matplotlib.pyplot as plt


def hdf5Analysis(filePaths,
                 attributeKey,
                 attributeName,
                 distinguishingFeature="Temperature",
                 title="",
                 single=False,
                 saveFig=False,
                 realTimeValue=False,
                 ylimLower=None,
                 ylimUpper=None):

    dataGroup = "/Dataset"
    stateGroup = dataGroup + "/State"
    settingsGroup = "/Setting"

    outputFiles         = []
    settings            = []
    distinguishers      = []
    desiredAttributes   = []
    timesteps           = []

    for i in range(0, len(filePaths)):
        outputFiles.append(h5py.File(filePaths[i]))
        settings.append(outputFiles[i][settingsGroup])

        if (distinguishingFeature == "Temperature"):
            distinguishers.append(settings[i].attrs["T"])
        elif(distinguishingFeature == "Particles"):
            volume = settings[i].attrs["nelx"] * settings[i].attrs["nely"] * settings[i].attrs["nelz"]
            distinguishers.append(volume)
        elif(distinguishingFeature == "volfrac"):
            distinguishers.append(settings[i].attrs["volfrac"])
        elif(distinguishingFeature == "timestep"):
            distinguishers.append(settings[i].attrs["dt"])
        else:
            raise Exception("Invalid Distinguishing Feature: " + distinguishingFeature + 
                            "/nValid features are: `Temperature` `Particles`")

        desiredAttributes.append(np.asarray(outputFiles[i][dataGroup + "/" + attributeKey]))
        timesteps.append(settings[i].attrs["dt"])


    print("Plotting", title)

    if (single):
        plotting.plotAttribute(timesteps[0],
                               desiredAttributes[0],
                               distinguishers[0],
                               attributeKey,
                               realTimeValue=realTimeValue,
                               ylimLower=ylimLower,
                               ylimUpper=ylimUpper,
                               savePlots=saveFig)
    else:
        plotting.plotAttributeSideBySide(desiredAttributes,
                                         distinguishers,
                                         distinguishingFeature,
                                         attributeName,
                                         timesteps,
                                         title=title,
                                         savePlots=saveFig,
                                         realTimeValue=realTimeValue,
                                         ylimLower=ylimLower,
                                         ylimUpper=ylimUpper)

    return

def hdf5AnalyzeSpecificAttributes(filePath,
                                  attributeKeys,
                                  attributeNames,
                                  title="",
                                  saveFig=False,
                                  xlims=None,
                                  ylims=None,
                                  annotation="",
                                  biggerSpace=False,
                                  poster=False,
                                  scale=1,
                                  png=False):

    dataGroup = "/Dataset"
    settingsGroup = "/Setting"

    distinguishers      = []
    desiredAttributes   = []
    
    outputFile  = h5py.File(filePath)
    settings    = outputFile[settingsGroup]

    for i in range(0, len(attributeKeys)):

        desiredAttributes.append(np.asarray(outputFile[dataGroup + "/" + attributeKeys[i]]))
        # timesteps.append(settings[i].attrs["dt"])


    print("Plotting", title)

    subtitles = []
    for name in attributeNames:
        subtitles.append(name + " at ")

    plotting.plotAttributesOnTop(   desiredAttributes,
                                    attributeNames,
                                    [title],
                                    ylims,
                                    xlims,
                                    savePlots=saveFig,
                                    markers=[True]*len(desiredAttributes),
                                    name=title,
                                    annotation=annotation,
                                    biggerSpace=biggerSpace,
                                    poster=poster,
                                    shareX=True,
                                    scale=scale,
                                    png=png)
                                    # linewidth=2)

    return

def overShootExample():
    filePath = ['../outputs/finally_good_for_real/10 Deg/hypopt_output_small_10deg_0.001dt_OVERSHOOT_EXAMPLE.h5']#,'../outputs/cyclone/hypopt_output_big_0.12vf_0.1deg_0.001ts_1000.h5']

    hdf5AnalyzeSpecificAttributes(filePath[0],
                                ["Temperature", "Compliance"],
                                ["Temperature", "Compliance"],
                                "Temperature 10 with Timestep 0.001",
                                ylims=[[-2, 20], [-8e-9, 40]],#1.1e11]],
                                xlims=[-5, 300],
                                saveFig=True,
                                # annotation="A",
                                biggerSpace=True)

def overShootExampleLowerTimestep():
    filePath = '../outputs/finally_good_for_real/10 Deg/hypopt_output_small_10deg_0.0001dt_OVERSHOOT_BETTER.h5'

    hdf5AnalyzeSpecificAttributes(filePath,
                                ["Temperature", "Compliance"],
                                ["Temperature", "Compliance"],
                                "Temperature 10 with Timestep 0.0001",
                                ylims=[[-2, 20], [-8e-9, 40]],#1.1e11]],
                                xlims=[-5, 3000],
                                saveFig=True,
                                # annotation="B",
                                biggerSpace=True)

def compareTwoFilesAllValues(filepaths):
    hdf5Analysis(filepaths, "Temperature",       "Temperature"           , distinguishingFeature="Temperature")#, single=True)
    hdf5Analysis(filepaths, "Lambda",            "Lagrange Multiplier" , distinguishingFeature="Temperature")#, single=True)
    hdf5Analysis(filepaths, "Hamiltonian",       "Hamiltonian"           , distinguishingFeature="Temperature")#, single=True)
    hdf5Analysis(filepaths, "Compliance",        "Compliance"            , distinguishingFeature="Temperature")#, single=True)
    hdf5Analysis(filepaths, "Volume Fraction",   "Volume Fraction"       , distinguishingFeature="Temperature")#, single=True)
    hdf5Analysis(filepaths, "Max Position",      "Max Position"          , distinguishingFeature="Temperature")#, single=True)

def analyzeDirectOutput(filePath=None, realTimeValue=False):
    if (None==filePath):
        filePath = ['../bin/hypopt_output.h5']
    else:
        filePath = [filePath]

    hdf5Analysis(filePath, "Temperature",       "Temperature"           , single=True, realTimeValue=realTimeValue)
    hdf5Analysis(filePath, "Lambda",            "Lagrange Multiplier" , single=True, realTimeValue=realTimeValue)
    hdf5Analysis(filePath, "Hamiltonian",       "Hamiltonian"           , single=True, realTimeValue=realTimeValue)
    hdf5Analysis(filePath, "Compliance",        "Compliance"            , single=True, realTimeValue=realTimeValue)
    hdf5Analysis(filePath, "Volume Fraction",   "Volume Fraction"       , single=True, realTimeValue=realTimeValue)
    hdf5Analysis(filePath, "Max Position",      "Max Position"          , single=True, realTimeValue=realTimeValue)

def multiFileSingleAttributeAnalysis(filePaths,
                                     attributeName,
                                     xlims=None,#[-10, 3000],
                                     ylims=None,
                                     logs=None,
                                     realTimeValue=False,
                                     markers=None,
                                     saveFig=False,
                                     desiredLines=None,
                                     linewidth=3.4):
    dataGroup = "/Dataset"
    settingsGroup = "/Setting"

    distinguishers      = []
    desiredAttributes   = []
    timesteps           = []
    temperatures        = []


    for i in range(0, len(filePaths)):
        outputFile  = h5py.File(filePaths[i])
        settings    = outputFile[settingsGroup]
        desiredAttributes.append(np.asarray(outputFile[dataGroup + "/" + attributeName]))
        timesteps.append(settings.attrs["dt"])
        temperatures.append(settings.attrs["T"])
        # timesteps.append(settings[i].attrs["dt"])


    print("Plotting ", attributeName, " analysis")

    subtitles = []
    for i in range(0, len(filePaths)):
        subtitles.append("Temperature " + str(temperatures[i]) + " with timestep " + str(timesteps[i]))

    plotting.plotAttributesOnTop(   desiredAttributes,
                                    [attributeName],#*len(filePaths),
                                    subtitles,
                                    ylims,
                                    xlims,
                                    savePlots=saveFig,
                                    realTimeValue=realTimeValue,
                                    timesteps=timesteps,
                                    logplots=logs,
                                    markers=markers,
                                    name=attributeName + "_anlaysis",
                                    desiredLines=desiredLines,
                                    linewidth=linewidth)


    return

# filePath2 = '../outputs/hypopt_output_first_run_broke.h5'
# filePath1 = '../outputs/hypopt_output_smallscale_system.h5'

# filePaths = [filePath1, filePath2]

# hdf5Analysis(filePaths, "Temperature", "Temperature", "Particles", "System Temperature at 0 Degrees", saveFig=True)

# filePath = ['../bin/hypopt_output.h5']#, '../bin/hypopt_output__0deg__0.01step.h5']
# filePath = ['../bin/hypopt_output__0deg__0.01step.h5']
# filePath = ['../outputs/finally_good_for_real/cyclone/hypopt_output_0.01deg_0.005dt_500.h5']
# filePath = ['../outputs/cyclone/hypopt_output_big_0.12vf_0.1deg_0.001ts_1000.h5']

# filePath = ['../outputs/finally_good_for_real/hypopt_output_10deg_0.001dt_OVERSHOOT_EXAMPLE.h5']#,'../outputs/cyclone/hypopt_output_big_0.12vf_0.1deg_0.001ts_1000.h5']

# hdf5AnalyzeSpecificAttributes(filePath[0],
#                             ["Temperature", "Compliance"],
#                             ["Temperature", "Compliance"],
#                             "5 Degree Simulation with 0.001 Timestep",
#                             ylims=[[-2, 20], [-8e-9, 40]],#1.1e11]],
#                             xlims=[-5, 300],
#                             saveFig=True)

# filePath = ['../outputs/finally_good_for_real/hypopt_output_small_0.001deg_0.01dt_50000.h5',
#             '../outputs/finally_good_for_real/hypopt_output_small_0deg_0.1dt_49000.h5']
# compareTwoFilesAllValues(filePath)
# filePath = "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/5 Deg/hypopt_output_small_5deg_0.00005dt_20000.h5"
# filePath = "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/hypopt_output_small_10deg_0.00001dt_25000.h5"
# # # analyzeDirectOutput(filePath, realTimeValue=False)
# analyzeDirectOutput(filePath, realTimeValue=True)
# overShootExampleLowerTimestep()


################## VALIDATION RESULTS ####################
if (False):

    analysisFiles = ["G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/0 Deg/hypopt_output_small_0deg_0.1dt_49000.h5",
                    "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/0.1 Deg/hypopt_output_small_0.1deg_0.0005deg_6000.h5",
                    "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/5 Deg/hypopt_output_small_5deg_0.00001dt_60000.h5"]#,
                    # "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/10 Deg/hypopt_output_small_10deg_0.00001dt_25000.h5"]
                    #  ?"G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/10 Deg/hypopt_output_small_10deg_0.0001dt_OVERSHOOT_BETTER.h5"]

    hamiltonianylims = [(5e-2, 10),
                        (-1, 600),
                        (-1, 33000),
                        (-1, 63000)]

    temperatureylims = [(5e-10, 0.001),
    # temperatureylims = [(-0.0001, 0.0005),
                        (-0.01, 0.15),
                        (-1, 10),
                        (-1, 20)]

    logsHamiltonian = [True, False, False, False]
    logsTemperature = logsHamiltonian
    markersHamiltonian = [False, False, False, False]
    markersTemperature = markersHamiltonian

    multiFileSingleAttributeAnalysis(analysisFiles,
                                    "Hamiltonian",
                                    realTimeValue=False,
                                    logs=logsHamiltonian,
                                    ylims=hamiltonianylims,
                                    saveFig=True,
                                    markers=markersHamiltonian)


    multiFileSingleAttributeAnalysis(analysisFiles,
                                    "Temperature",
                                    realTimeValue=False,
                                    logs=logsTemperature,
                                    ylims=temperatureylims,
                                    saveFig=True,
                                    markers=markersTemperature,
                                    desiredLines=[0, 0.1, 5, 10])

    multiFileSingleAttributeAnalysis(analysisFiles,
                                    "Volume Fraction",
                                    realTimeValue=False,
                                    #  logs=logsTemperature,
                                    #  ylims=[(0.1199999,0.120001)]*4)
                                    saveFig=True,
                                    linewidth=1.5)
                                    #  markers=markersTemperature,
                                    #  desiredLines=[0.12, 0, 0, 0])

    multiFileSingleAttributeAnalysis(analysisFiles,
                                 "Compliance",
                                 realTimeValue=False,
                                 saveFig=True
                                 )

################## Compliance Oscillations ####################

if (False):
    filePath = ["../outputs/finally_good_for_real/0.001 Deg/hypopt_output_small_0.001deg_0.01dt_50000.h5"]

    hdf5Analysis(filePath,
                "Compliance",
                "Compliance",
                single=True,
                realTimeValue=False,
                ylimLower=0,
                ylimUpper=4,
                saveFig=True)

################## Overshoot ####################
if (False):
    overShootExample()
    overShootExampleLowerTimestep()


################## Poster Analasys ####################
if(False):
    filePath = "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/0.1 Deg/hypopt_output_small_0.1deg_0.0005deg_6000.h5"

    attributes  = ["Temperature", "Hamiltonian", "Volume Fraction"]
    names       = ["Temp", "Hamilt", "Vol Frac"]


    hdf5AnalyzeSpecificAttributes(filePath,
                                attributes,
                                names,
                                "Desired Temperature 0.1",
                                ylims=[[-0.02,0.17], [-50, 700], [0.12-0.5e-7, 0.12+6e-7]],#1.1e11]],
                                # xlims=[-5, 300],
                                saveFig=True,
                                # annotation="A",
                                biggerSpace=True,
                                poster=True,
                                scale=0.7,
                                png=True)

if (True):
    filePath = "../outputs/finally_good_for_real/0.001 Deg/hypopt_output_small_0.001deg_0.01dt_50000.h5"

    attributes  = ["Compliance"]
    names       = ["Compliance"]


    hdf5AnalyzeSpecificAttributes(filePath,
                                attributes,
                                names,
                                "Desired Temperature 0.001",
                                ylims=[[0,2.5]],
                                xlims=[-50, 30000],
                                saveFig=True,
                                # annotation="A",
                                # biggerSpace=True,
                                poster=True,
                                scale=1,
                                png=True)

################## Specific File ####################
if(False):
#    filePath = "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/5 Deg/hypopt_output_small_5deg_0.00001dt_60000.h5"
    filePath = "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/cyclone/hypopt_output_big_0deg_0.1dt_12000_restart3_lastpls.h5"
    # filePath = "G:/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/0.001 Deg/hypopt_output_small_0.001deg_0.01dt_50000.h5"
    analyzeDirectOutput(filePath, realTimeValue=True)


################## ./bin ####################
if(False):
    analyzeDirectOutput(realTimeValue=False)

plt.show()