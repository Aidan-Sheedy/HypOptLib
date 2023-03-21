import matplotlib.pyplot as plt
import numpy as np

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
        save_name = "../figures/" + str(ylabel) + ".pdf"
        figure.savefig(save_name, format='pdf', dpi=1200,bbox_inches = 'tight')

    return