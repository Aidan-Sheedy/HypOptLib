import matplotlib.pyplot as plt
import numpy as np
import os

import config

def plotAttribute(timestep, attribute, temperature, ylabel, ylimLower=0, ylimUpper=None, savePlots=False):
    '''
    Plots a given attribute against time given the timestep.
    '''

    numPoints = len(attribute)
    time = np.arange(0,numPoints)# * timestep

    figure = plt.figure()
    figure.set_size_inches(5, 4)

    figure.suptitle(ylabel + " at Temperature " + str(temperature))

    plt.plot(time, attribute, marker='.')
    plt.xlabel("timestep")
    plt.ylabel(ylabel)

    if (ylimUpper != None):
        plt.ylim(ylimLower, ylimUpper)


    plt.grid(linewidth=1.3, linestyle=':')

    # plt.show()

    if (savePlots):
        save_name = "../figures/" + ylabel + " at Temperature " + str(temperature) + ".pdf"
        figure.savefig(save_name, format='pdf', dpi=1200,bbox_inches = 'tight')

    return


def plotAttributeSideBySide(attributes, distinguisingFeatures, distingusingFeatureName, ylabel, title="", timestep=0.001, ylimLower=0, ylimUpper=None, savePlots=False, units=""):
    '''
    Plots a given attribute against time given the timestep.
    '''
    numPlots = len(attributes)

    time = []*numPlots

    for i in range(0, numPlots):
        time.append(np.arange(0,len(attributes[i])))# * timestep

    figure, subplots = plt.subplots(1, numPlots)
    figure.set_size_inches(5*numPlots, 4)

    if ("" != title):
        figure.suptitle(title, fontsize=config.figureTitleFontsize)

    colours = ['b', 'g', 'r', 'tab:orange']

    for i in range(0,numPlots):
        subplots[i].set_title(str(distinguisingFeatures[i]) + " " + distingusingFeatureName)#, fontsize=config.generalFontsize)
        subplots[i].plot(time[i], attributes[i], marker='.', c=colours.pop(0))
        subplots[i].set_xlabel("Timestep")#, fontsize=config.generalFontsize)
        if (i == 0):
            if ("" != units):
                subplots[i].set_ylabel(ylabel + " [" + units + "]")#, fontsize=config.generalFontsize)
            else:
                subplots[i].set_ylabel(ylabel)#, fontsize=config.generalFontsize)

        subplots[i].grid(linewidth=1.3, linestyle=':')
        subplots[i].tick_params(axis='both', labelsize=config.axisTickFontsize)

    figure.tight_layout(pad=0.3)
    figure.subplots_adjust(top=0.8)

    # plt.show()

    if (savePlots):
        if (not os.path.isdir("../figures")):
            os.mkdir("../figures")
        save_name = "../figures/" + ylabel + "_comparisons_at_" + distingusingFeatureName
        for feature in distinguisingFeatures:
            save_name += "_" + str(feature)
        save_name += ".pdf"
        figure.savefig(save_name, format='pdf', dpi=1200,bbox_inches = 'tight')

    return