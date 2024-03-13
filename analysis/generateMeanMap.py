#!/usr/bin/env python3

'''
This file calculates the mean of positions accross all converged iterations. It can support
an arbitrary number of restart files.

Author: Aidan Sheedy

Copyright (C) 2024 Aidan Sheedy

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
'''

import h5py
import os
import numpy as np
from utilities import *

files = ['SampleFiles/sample_cantilever.h5',
         'SampleFiles/sample_cantilever_restart.h5']

# Manually decide when the simulation has converged.
converged_iteration = 8580
useIterationConverged = False
converged_time = 0.6

# Get settings from file. If multiple files, assume each file has the same settings as the first.
hdf5File = h5py.File(files[0])
nelx = hdf5File["/Setting"].attrs["nelx"]
nely = hdf5File["/Setting"].attrs["nely"]
nelz = hdf5File["/Setting"].attrs["nelz"]
dt = hdf5File["/Setting"].attrs["dt"]
hdf5File.close()

if useIterationConverged:
    converged_time = converged_iteration*dt

# Initialize variables
mean_map = np.zeros(shape=(nelz, nely, nelx))
num_iterations = 0
i = 0

# Sum position values in each file
for file in files:
    hdf5File = h5py.File(file)

    try:
        saveFreq = hdf5File["/Setting"].attrs["saveFreq"]
    except:
        print("Provided file does not indicate save frequency. Using 10 as default.")
        saveFreq = 10

    names = getSortedNames(hdf5File["/Dataset/State"])
    for name in names:
        # The iterations start saving at iteration 0.
        # So iteration 0 will have timestamp = dt, but the next timestamp will be
        # dt * (saveFreq + 1), and so on.
        timestamp = dt * (i*saveFreq + 1)
        i += 1

        # only sum iterations that have converged.
        if ( converged_time < timestamp ):
            iteration = np.asarray(hdf5File[name])
            mean_map += np.asarray(iteration)
            num_iterations += 1
    hdf5File.close()

# Divide by number of iterations.
mean_map = mean_map / num_iterations

# Create an hdf5 and xmf file with the mean map information.
new_file_basename = os.path.splitext(os.path.basename(files[0]))[0] + "_mean_map"
new_filename =  os.path.dirname(files[0]) + "/" + new_file_basename
print("Generating mean map to: ", new_filename + ".xmf")

# hdf5 file:
map_file = h5py.File(new_filename + ".h5", 'w')
map_file.create_dataset("Position_Mean_Map", data=mean_map)
map_file.close()

# xmf file
xdmfFile = xmf(new_filename + ".xmf")
xdmfFile.startGenericFile()
xdmfFile.insertRectilinearMesh([nelx+1, nely+1, nelz+1])
xdmfFile.insertAttribute(new_file_basename + ".h5", "Mean Position", "Position_Mean_Map")
xdmfFile.closeGenericFile()
