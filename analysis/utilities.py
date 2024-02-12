#!/usr/bin/env python3

'''
This utilities file provides generic utilities useful for analyzing HypOptLib output files.

Author: Aidan Sheedy

TODO This file needs license information.
'''


def getSortedNames(iterationNames):
    '''
    Sorts a set of HypOptLib iteration names from lowest to highest.

    Parameters
    ----------
    iterationNames : list of names to sort

    Returns
    -------
    Sorted list of names, with 'Dataset/State/' pre-pended to each iteration.
    '''

    i = 0
    sorted = [0]*len(iterationNames)
    for name in iterationNames:
        sorted[i] = int(name[9:])
        i += 1
    sorted.sort()
    i = 0
    for name in sorted:
        sorted[i] = 'Dataset/State/iteration' + str(name)
        i+=1

    return sorted

class xmf:
    '''
    Class designed to handle creating HypOptLib specific xdmf files to interface with ParaView.

    Methods
    -------
    startGenericFile()
        Writes the header for an arbitrary xmf file. Must be closed with closeGenericFile().

    closeGenericFile()
        Writes the footer to a generic xmf file.

    insertRectilinearMesh(dimensions)
        Inserts a rectilinear mesh of the provided dimensions.

    insertAttribute(fileName, attributeName, attributeLocation)
        Inserts the desired attribute into the xmf file.

    startTimesteppedFile()
        Writes the header for a timestepped xmf file. Must be closed with closeTimesteppedFile().

    closeTimesteppedFile()
        Writes the footer to a timestepped xmf file.

    insertTimestep(time, hdf5FileName, iterationName)
        Inserts a timestep into an opened timestepped xmf file.
    '''

    def __init__(self, filepath):
        self.__filepath = filepath

        # Write header
        with open(self.__filepath, 'w') as file:
            file.write("<?xml version=\"1.0\" ?>\n")
            file.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
            file.write("<Xdmf Version=\"2.0\">\n")

    def startGenericFile(self):
        '''
        Writes the header for an arbitrary xmf file. Must be closed with closeGenericFile().
        '''
        with open(self.__filepath, 'a') as file:
            file.write("\t<Domain>\n")
            file.write("\t\t<Grid Name =\"mesh\" GridType=\"Uniform\">\n")

    def insertRectilinearMesh(self, dimensions):
        '''
        Inserts a rectilinear mesh of the provided dimensions.

        Parameters
        ----------
        dimensions : the (x, y, z) dimensions of the mesh.
        '''
        self.__dimensions = dimensions
        with open(self.__filepath, 'a') as file:
            file.write("\t\t\t<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" + str(dimensions[2]) + " " + str(dimensions[1]) + " "+ str(dimensions[0]) + "\"/>\n")
            file.write("\t\t\t<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n")
            file.write("\t\t\t\t<DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">0.0 0.0 0.0</DataItem>\n")
            file.write("\t\t\t\t<DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">1.0 1.0 1.0</DataItem>\n")
            file.write("\t\t\t</Geometry>\n")

    def insertAttribute(self, hdf5FileName, attributeName, attributeLocation):
        '''
        Inserts the desired attribute into the xmf file.

        Parameters
        ----------
        hdf5FileName : the name of the hdf5 file in which the timestep is stored. 

        attributeName : the name of the attribute

        attributeLocation : the location of the attribute in the hdf5 file.
        '''
        with open(self.__filepath, 'a') as file:
            file.write("\t\t\t<Attribute Name=\"" + attributeName + "\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
            file.write("\t\t\t\t<DataItem Dimensions=\"" + str(self.__dimensions[2]-1) + " " + str(self.__dimensions[1]-1) + " " + str(self.__dimensions[0]-1) +
                        "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n")
            file.write("\t\t\t\t " + hdf5FileName + ":" + attributeLocation + "\n")
            file.write("\t\t\t\t</DataItem>\n")
            file.write("\t\t\t</Attribute>\n")

    def closeGenericFile(self):
        '''
        Writes the footer to a generic xmf file.
        '''
        with open(self.__filepath, 'a') as file:
            file.write("\t\t</Grid>\n")
            file.write("\t</Domain>\n")
            file.write("</Xdmf>\n")

    def startTimesteppedFile(self, dimensions):
        '''
        Writes the header for a timestepped xmf file. Must be closed with closeTimesteppedFile().
        '''
        self.__dimensions = dimensions

        with open(self.__filepath, 'a') as file:
            file.write("\t<Domain>\n")
            file.write("\t\t<Grid Name =\"mesh\" GridType=\"Collection\" CollectionType=\"Temporal\">\n\n")

    def insertTimestep(self, time, hdf5FileName, iterationName):
        '''
        Inserts a timestep into an opened timestepped xmf file.

        The timestep can be any arbitrary value, and does not need to be consecutive with the previous timestep.
        However, it should be in order from smallest to largest.

        Parameters
        ----------
        time : the timestamp for the timestep.

        hdf5FileName : the name of the hdf5 file in which the timestep is stored.

        iterationName : the location of the iteration in the hdf5 file.
        '''
        with open(self.__filepath, 'a') as file:
            file.write("\t\t<Grid Name =\"mesh\" Type=\"Uniform\">\n")
            file.write("\t\t<Time Type =\"Single\" Value=\""+str(time)+"\"/>\n")
            file.write("\t\t\t<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\""+str(self.__dimensions[2]+1)+" "+str(self.__dimensions[1]+1)+" "+str(self.__dimensions[0]+1)+"\"/>\n")
            file.write("\t\t\t<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n")
            file.write("\t\t\t\t<DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">0.0 0.0 0.0</DataItem>\n")
            file.write("\t\t\t\t<DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">1.0 1.0 1.0</DataItem>\n")
            file.write("\t\t\t</Geometry>\n")
            file.write("\t\t\t<Attribute Name=\"Position\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
            file.write("\t\t\t\t<DataItem Dimensions=\""+str(self.__dimensions[2])+" "+str(self.__dimensions[1])+" "+str(self.__dimensions[0])+"\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n")
            file.write("\t\t\t\t" + hdf5FileName + ":" + iterationName + "\n")
            file.write("\t\t\t\t</DataItem>\n")
            file.write("\t\t\t</Attribute>\n")
            file.write("\t\t\t</Grid>\n\n")

    def closeTimesteppedFile(self):
        '''
        Writes the footer to a timestepped xmf file.
        '''
        with open(self.__filepath, 'a') as file:
            file.write("\t\t</Grid>\n")
            file.write("\t</Domain>\n")
            file.write("</Xdmf>\n")

