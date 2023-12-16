#!/usr/bin/env python3

import h5py
import numpy as np

fileName = "hypopt_output_big_0deg_0.01dt_4319"

# Copy the desired data set 
filePath = '/mnt/g/Projects/ENPH455/PETSc_paper_github/Hyperoptimization_using_Petsc/outputs/finally_good_for_real/big_0deg/' + fileName + '.h5'
dataGroup = '/Dataset'


ogFile = h5py.File(filePath)
tempFile = h5py.File("../xmfFiles/" + fileName + ".h5", "w")

tempFile["Final Position"] = np.asarray(ogFile[dataGroup + "/State/iteration4319"])

ogFile.close()
tempFile.close

# Create the xdmf file
xdmfOutputPath = "../xmfFiles/" + fileName + ".xmf"

xdmfFile = open(xdmfOutputPath, "w")

# Write the header
xdmfFile.write("<?xml version=\"1.0\" ?>\n")
xdmfFile.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
xdmfFile.write("<Xdmf Version=\"2.0\">\n")

# Write the mesh description
xdmfFile.write("\t<Domain>\n")

xdmfFile.write("\t\t<Grid Name =\"mesh\" GridType=\"Uniform\">\n")
xdmfFile.write("\t\t\t<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"65 65 128\"/>\n")
xdmfFile.write("\t\t\t<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n")
xdmfFile.write("\t\t\t\t<DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">0.0 0.0 0.0</DataItem>\n")
xdmfFile.write("\t\t\t\t<DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">1.0 1.0 1.0</DataItem>\n")
xdmfFile.write("\t\t\t</Geometry>\n")
xdmfFile.write("\t\t\t<Attribute Name=\"Position\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
xdmfFile.write("\t\t\t\t<DataItem Dimensions=\"16 16 32\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n")
xdmfFile.write("\t\t\t\t " + fileName + ".h5:/Final Position\n")
xdmfFile.write("\t\t\t\t</DataItem>\n")
xdmfFile.write("\t\t\t</Attribute>\n")
xdmfFile.write("\t\t</Grid>\n")
xdmfFile.write("\t</Domain>\n")
xdmfFile.write("</Xdmf>\n")
