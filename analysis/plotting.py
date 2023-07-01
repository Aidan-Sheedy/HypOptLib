import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os

import config

def plotAttribute(timestep, attribute, temperature, ylabel, ylimLower=0, ylimUpper=None, savePlots=False, realTimeValue=False):
    '''
    Plots a given attribute against time given the timestep.
    '''

    numPoints = len(attribute)
    time = np.arange(0,numPoints)
    if (realTimeValue):
        time = time * timestep

    figure = plt.figure()
    figure.set_size_inches(8, 2.5)

    figure.suptitle(ylabel + " at Temperature " + str(temperature))

    plt.plot(time, attribute, c='b')
    if (realTimeValue):
        plt.xlabel("Time")
    else:
        plt.xlabel("Iteration")    

    plt.ylabel(ylabel)

    if (ylimUpper != None):
        plt.ylim(ylimLower, ylimUpper)


    plt.grid(linewidth=1.3, linestyle=':')
    
    figure.subplots_adjust(top=0.85)

    # plt.show()

    if (savePlots):
        save_name = "../figures/" + ylabel + " at Temperature " + str(temperature) + ".pdf"
        figure.savefig(save_name, format='pdf', dpi=1200,bbox_inches = 'tight')

    return


def plotAttributeSideBySide(attributes, distinguisingFeatures, distingusingFeatureName, ylabel, timesteps, title="", ylimLower=0, ylimUpper=None, savePlots=False, units="", realTimeValue=False):
    '''
    Plots a given attribute against time given the timestep.
    '''
    numPlots = len(attributes)

    time = []*numPlots

    for i in range(0, numPlots):
        time.append(np.arange(0,len(attributes[i])))# * timestep
        if (realTimeValue):
            time[i] *= timesteps[i]

    figure, subplots = plt.subplots(1, numPlots)
    figure.set_size_inches(5*numPlots, 4)

    if ("" != title):
        figure.suptitle(title, fontsize=config.figureTitleFontsize)

    colours = ['b', 'g', 'r', 'tab:orange']

    for i in range(0,numPlots):
        subplots[i].set_title(str(distinguisingFeatures[i]) + " " + distingusingFeatureName)#, fontsize=config.generalFontsize)
        # subplots[i].plot(time[i], attributes[i], marker='.', c=colours.pop(0))
        subplots[i].plot(time[i], attributes[i], c=colours.pop(0))
        if (realTimeValue):
            subplots[i].set_xlabel("Time")
        else:
            subplots[i].set_xlabel("Iteration")#, fontsize=config.generalFontsize)
        if (i == 0):
            if ("" != units):
                subplots[i].set_ylabel(ylabel + " [" + units + "]")#, fontsize=config.generalFontsize)
            else:
                subplots[i].set_ylabel(ylabel)#, fontsize=config.generalFontsize)

        subplots[i].grid(linewidth=1.3, linestyle=':')
        subplots[i].tick_params(axis='both', labelsize=config.axisTickFontsize)

    figure.tight_layout(pad=0.4)
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

def plotAttributesOnTop(attributes,
                        ylabels,
                        titles,
                        ylims,
                        xlims,
                        savePlots=False,
                        units="",
                        markers=None,
                        timesteps=None,
                        realTimeValue=False,
                        logplots=None,
                        name="",
                        desiredLines=None,
                        linewidth=4,
                        annotation="",
                        biggerSpace=False,
                        biggestSpace=True,
                        poster=False,
                        shareX=False,
                        scale=1,
                        png=False):
    '''
    Plots a given attribute against time given the timestep.
    '''

    numPlots = len(attributes)
    print("######\nNUMPLOTS: ", numPlots)

    if (None == logplots):
        logplots = [False]*numPlots

    if (None == markers):
        markers = [False]*numPlots

    time = []*numPlots

    for i in range(0, numPlots):
        time.append(np.arange(0,len(attributes[i])))# * timesteps[i]
        if (realTimeValue):
            time[i] = time[i] * timesteps[i]

    if (shareX):
        figure, subplots = plt.subplots(numPlots, sharex=True)
    else:
        figure, subplots = plt.subplots(numPlots)

    if (1 == numPlots):
        subplots = [subplots]
    if (biggerSpace):
        figure.set_size_inches(7*scale, 2*numPlots*scale)
    elif(3 == numPlots or biggestSpace):
        figure.set_size_inches(7*scale, 3*numPlots*scale)
    else:
        figure.set_size_inches(7*scale, 1.8*numPlots*scale)
    
    # figure.suptitle(title, fontsize=config.figureTitleFontsize)

    if (1 == len(ylabels)):
        figure.supylabel(ylabels[0])

    if (poster):
        colours = ['#e74d4a', '#b13b39', '#740F0E', 'r'] # RED
        # colours = ['#92b9f5', '#567ebb', '#2c3d55', 'r'] # BLUE
    else:
        colours = ['b', 'g', 'r', 'tab:orange']

    for i in range(0,numPlots):
        if ("" != annotation and i == 0):
            subplots[i].annotate(   annotation,
                                    xy=(0, 1),
                                    xycoords='figure fraction',
                                    horizontalalignment='left',
                                    verticalalignment='top',
                                    fontsize=config.figureTitleFontsize+10,
                                    color='#00007f')
        # subplots[i].set_title(titles[i])#, fontsize=config.generalFontsize)
        if (markers[i]):
            if (logplots[i]):
                subplots[i].semilogy(time[i], attributes[i], marker='.', c=colours.pop(0), markersize=4.5)
            else:
                subplots[i].plot(time[i], attributes[i], marker='.', c=colours.pop(0), markersize=4.5)
        else:
            if (logplots[i]):
                subplots[i].semilogy(time[i], attributes[i], c=colours.pop(0), linewidth=linewidth)
                if (np.abs((np.log10(ylims[i][1]) - np.log10(ylims[i][0]))) > 3 ):
                    subplots[i].set_yticks([10e-6, 10e-8, 10e-10])
                    # subplots[i].yaxis.set_major_locator(mticker.LogLocator(numticks=5, subs="auto"))
            else:
                subplots[i].plot(time[i], attributes[i], c=colours.pop(0), linewidth=linewidth)
        if (i == numPlots-1):
            if (realTimeValue):
                subplots[i].set_xlabel("Time")
            else:
                subplots[i].set_xlabel("Iteration")#, fontsize=config.generalFontsize)
        if (None != desiredLines):
            subplots[i].plot(time[i], [desiredLines[i]]*len(time[i]), c='k')
        if (1 < len(ylabels)):
            subplots[i].set_ylabel(ylabels[i])#, fontsize=config.generalFontsize)
        print("I!!!!! ", i)
        if (len(titles) > 1):
            subplots[i].set_title(titles[i], fontsize=config.generalFontsize)
        elif(len(titles) == 1 and i == 0):
            subplots[i].set_title(titles[0], fontsize=config.generalFontsize)

        # if (logplots[i] and (np.abs((np.log10(ylims[i][1]) - np.log10(ylims[i][0]))) > 3 ) ):
        #     subplots[i].yaxis.set_major_locator(mticker.MaxNLocator(4))

        subplots[i].grid(linewidth=1.3, linestyle=':')
        subplots[i].tick_params(axis='both', labelsize=config.axisTickFontsize)

        if (None != ylims):
            subplots[i].set_ylim(ylims[i][0], ylims[i][1])
        if (None != xlims):
            subplots[i].set_xlim(xlims[0], xlims[1])

    # figure.subplots_adjust(top=0.9)
    if (biggerSpace):
        figure.tight_layout(pad=0.6)
    else:
        figure.tight_layout(pad=0.6)
    figure.align_ylabels()

    # plt.show()

    if (savePlots):
        format="pdf"
        if(png):
            format="png"
        if (not os.path.isdir("../figures")):
            os.mkdir("../figures")
        save_name = "../figures/" + name + "." + format
        figure.savefig(save_name, format=format, dpi=1200,bbox_inches = 'tight')
        print("Saved figure:\n", save_name)

    return
