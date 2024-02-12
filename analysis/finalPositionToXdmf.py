#!/usr/bin/env python3

'''
This file provides a function to generate an xmf file which interfaces the final position of a
HypOptLib output file with ParaView.

Author: Aidan Sheedy

TODO This file needs license information.
'''

import h5py
import os
from utilities import *

# Copy the desired data set 
file = 'SampleFiles/sample_cantilever_restart.h5'

# Get the dimensions of the file
hdf5File = h5py.File(file)
mesh_dimensions = [hdf5File["/Setting"].attrs["nelx"]+1, hdf5File["/Setting"].attrs["nely"]+1, hdf5File["/Setting"].attrs["nelz"]+1]
hdf5File.close()

# The xdmf file should be in the same directory as the original hdf5 file.
filePath = os.path.dirname(file)
file_justname = os.path.splitext(os.path.basename(file))[0]
xdmfOutputPath = filePath + "/" + file_justname + "_final_position.xmf"

print("Writing xmf file:", xdmfOutputPath)

# Create the xdmf file
xmfFile = xmf(xdmfOutputPath)
xmfFile.startGenericFile()
xmfFile.insertRectilinearMesh(mesh_dimensions)
xmfFile.insertAttribute(file_justname + ".h5", "Final Position", "Dataset/Final Position")
xmfFile.closeGenericFile()
