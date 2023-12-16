#!/usr/bin/env python3

import h5py
import numpy as np

base = '../outputs/multicore_comparisons/32x16x16/multicore_deviation_6_mpi.h5'

dataGroup = '/Dataset'
baseFile = h5py.File(base)
basePositions = np.asarray(baseFile[dataGroup + "/" + "Final Position"])

average_min = 0
average_max = 0
average_mean = 0

print("Min point: ", np.min(basePositions))

for i in range(2, 11):
    compare = '../run/multicore_deviation_6_mpi_' + str(i) + '.h5'

    compareFile = h5py.File(compare)
    comparePositions = np.asarray(compareFile[dataGroup + "/" + "Final Position"])

    percentDifference = np.divide((basePositions - comparePositions), basePositions) * 100
 
    for i in range(0, percentDifference.shape[0]):
        for j in range(0, percentDifference.shape[1]):
            for k in range(0, percentDifference.shape[2]):
                if (percentDifference[i][j][k] < 0):
                    percentDifference[i][j][k] = -percentDifference[i][j][k]

    min = np.min(percentDifference)
    max = np.max(percentDifference)
    mean = np.mean(percentDifference)
    print("\nMin point: ", np.min(percentDifference))
    print("Percent difference min: ", min, "\nPercent Difference max: ", max, "\nPercent difference mean: ", mean)


    average_min = average_min + min
    average_max = average_max + max
    average_mean = average_mean + mean

average_min = average_min / 9
average_max = average_max / 9
average_mean = average_mean / 9


print("\nAverage Min: ", average_min, "\nAverage Max: ", average_max, "\nAverage Mean: ", average_mean)