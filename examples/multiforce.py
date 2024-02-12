#!/usr/bin/env python3

'''
This file runs a multi-force beam hyperoptimization simulation.

Author: Aidan Sheedy

TODO This file needs license information.
'''

import HypOptLib

######################################################################
# Initialize variables
solver = HypOptLib.HypOptLib()
domain = HypOptLib.DomainCoordinates()
fixedPoints     = HypOptLib.BoundaryCondition()
zforceCentre    = HypOptLib.BoundaryCondition()
zforceCorner1   = HypOptLib.BoundaryCondition()
zforceCorner2   = HypOptLib.BoundaryCondition()
yforceCentre    = HypOptLib.BoundaryCondition()
yforceCorner1   = HypOptLib.BoundaryCondition()
yforceCorner2   = HypOptLib.BoundaryCondition()

######################################################################
# Setup domain. This describes a rectangular prism of dimensions 2x1x1, with cubic voxels.
domain.xMinimum = 0
domain.xMaximum = 2
domain.yMinimum = 0
domain.yMaximum = 1
domain.zMinimum = 0
domain.zMaximum = 1

solver.setGridProperties([32, 16, 16], domain)

######################################################################
# Set up boundary conditions

# First boundary condition fixes the x=0 plane
fixedPoints.degreesOfFreedom = {0, 1, 2}
fixedPoints.type    = HypOptLib.BoundaryConditionType.FIXED_POINT
fixedPoints.xRange  = [0, 0]
fixedPoints.yRange  = [0, 1]
fixedPoints.zRange  = [0, 1]
fixedPoints.value   = 0

#---------------------------------------------------------------------
# Z Load 

# Set a line force at X=1, Z=0 in the Z DOF
zforceCentre.degreesOfFreedom = {2}
zforceCentre.type    = HypOptLib.BoundaryConditionType.LOAD
zforceCentre.xRange  = [2, 2]
zforceCentre.yRange  = [0, 1]
zforceCentre.zRange  = [0, 0]
zforceCentre.value   = -0.001

# Set the (1,0,0) corner to be half the line force
zforceCorner1.degreesOfFreedom = {2}
zforceCorner1.type    = HypOptLib.BoundaryConditionType.LOAD
zforceCorner1.xRange  = [2, 2]
zforceCorner1.yRange  = [0, 0]
zforceCorner1.zRange  = [0, 0]
zforceCorner1.value   = -0.0005

# Set the (1,1,0) corner to be half the line force
zforceCorner2.type    = HypOptLib.BoundaryConditionType.LOAD
zforceCorner2.degreesOfFreedom = {2}
zforceCorner2.xRange  = [2, 2]
zforceCorner2.yRange  = [1, 1]
zforceCorner2.zRange  = [0, 0]
zforceCorner2.value   = -0.0005

#---------------------------------------------------------------------
# y Load 

# Set a line force at X=1, y=1 in the Z DOF
zforceCentre.degreesOfFreedom = {2}
zforceCentre.type    = HypOptLib.BoundaryConditionType.LOAD
zforceCentre.xRange  = [2, 2]
zforceCentre.yRange  = [1, 1]
zforceCentre.zRange  = [0, 1]
zforceCentre.value   = 0.001

# Set the (1,1,1) corner to be half the line force
zforceCorner1.degreesOfFreedom = {2}
zforceCorner1.type    = HypOptLib.BoundaryConditionType.LOAD
zforceCorner1.xRange  = [2, 2]
zforceCorner1.yRange  = [1, 1]
zforceCorner1.zRange  = [1, 1]
zforceCorner1.value   = 0.0005

# Sets the (1,1,0) corner to be half the line force
zforceCorner2.type    = HypOptLib.BoundaryConditionType.LOAD
zforceCorner2.degreesOfFreedom = {2}
zforceCorner2.xRange  = [2, 2]
zforceCorner2.yRange  = [1, 1]
zforceCorner2.zRange  = [0, 0]
zforceCorner2.value   = 0.0005

# Apply boundary conditions
solver.setBoundaryConditions([fixedPoints, zforceCentre, zforceCorner1, zforceCorner2,
                              yforceCentre, yforceCorner1, yforceCorner2])

######################################################################
# Setup Hyperoptimization parameters.
solver.setSavePath("multiforce_beam.h5")
solver.setTargetTemperature(0.01)
solver.setTimestep(0.01)
solver.setMaximumIterations(6000)

######################################################################
# Run the simulation
saveRange = [0,6000]
solver.newRun(saveRange)
