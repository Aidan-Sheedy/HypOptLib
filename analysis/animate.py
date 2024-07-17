#!/usr/bin/env python3

'''
This file provides a function to generate an xmf file which interfaces a
HypOptLib output file with ParaView.

The function supports passing multiple restarted files and will order the timesteps
by the order of the files provided. The function assumes constant timestepped files.

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
from utilities import *

def generateTimestepXmf(files, saveFreq=0):
    # Parse if one or multiple files are used
    if (type(files) is list) or (type(files) is tuple):
        print("Generating XMF file for multiple files: ", files)
    else:
        print("Generating XMF file for one file: ", files)
        files = [files]

    # For each file, get the list of iteration names and other metadata needed for the xmf file.
    fileList = []
    firstFile = True
    for file in files:
        file_justname = os.path.splitext(os.path.basename(file))[0]
        print("----- parsing file: ", file)

        # Get ordered list of iterations
        original_file = h5py.File(file)
        dataset = original_file["/Dataset/State"]
        iterationNames = list(dataset.keys())
        iterationNames = getSortedNames(iterationNames)

        # Assume timestep is the same for each file.
        timestep = original_file["/Setting"].attrs["dt"]

        # Earlier files did not save the frequency option. If this is passed in, don't
        # try to parse it from the hdf5 file.
        if 0 == saveFreq:
            saveFreq = timestep = original_file["/Setting"].attrs["saveFreq"]

        # Make sure that each file has matching mesh
        if firstFile:
            nelx = original_file["/Setting"].attrs["nelx"]
            nely = original_file["/Setting"].attrs["nely"]
            nelz = original_file["/Setting"].attrs["nelz"]
            firstFile = False
        else:
            if ( (nelx != original_file["/Setting"].attrs["nelx"]) or
                 (nely != original_file["/Setting"].attrs["nely"]) or
                 (nelz != original_file["/Setting"].attrs["nelz"]) ):
                print("ERROR! The dimensions of each file provided do not match.")
                original_file.close()
                return

        # Setup list of iterations along with file name
        for name in iterationNames:
            fileList.append((file_justname + ".h5", name))

        original_file.close()

    # The xmf file must be created in the same directory as the original hdf5 files
    file_directory = os.path.dirname(files[0])
    filename = os.path.splitext(os.path.basename(files[0]))[0]
    new_path = str(file_directory)
    xdmfFileName  = new_path + "/" + filename + ".xmf"

    # Generate xmf file
    xdmfFile = xmf(xdmfFileName)
    xdmfFile.startTimesteppedFile([nelx, nely, nelz])

    # Iterate over each timestep
    i = 0
    for iteration in fileList:
        time = round(timestep * (i*saveFreq + 1), len(str(timestep)))
        xdmfFile.insertTimestep(time, iteration[0], iteration[1])
        i += 1

    # Write the footer
    xdmfFile.closeTimesteppedFile()

files = ['SampleFiles/sample_cantilever.h5',
         'SampleFiles/sample_cantilever_restart.h5']

generateTimestepXmf(files, 10)

