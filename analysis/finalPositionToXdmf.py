#!/usr/bin/env python3

'''
This file provides a function to generate an xmf file which interfaces the final position of a
HypOptLib output file with ParaView.

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
