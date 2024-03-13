#!/usr/bin/env python3

'''
This file runs a basic MBB beam hyperoptimization simulation.

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

import HypOptLib

######################################################################
# Initialize variables
solver = HypOptLib.HypOptLib()
domain = HypOptLib.DomainCoordinates()
fixedPoints = HypOptLib.BoundaryCondition()
forceCentre = HypOptLib.BoundaryCondition()
forceCorner1  = HypOptLib.BoundaryCondition()
forceCorner2  = HypOptLib.BoundaryCondition()

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

# Second boundary condition sets a line force at X=1, Z=0 in the Z DOF
forceCentre.degreesOfFreedom = {2}
forceCentre.type    = HypOptLib.BoundaryConditionType.LOAD
forceCentre.xRange  = [2, 2]
forceCentre.yRange  = [0, 1]
forceCentre.zRange  = [0, 0]
forceCentre.value   = -0.001

# Third boundary condition sets the (1,0,0) corner to be half the line force
forceCorner1.degreesOfFreedom = {2}
forceCorner1.type    = HypOptLib.BoundaryConditionType.LOAD
forceCorner1.xRange  = [2, 2]
forceCorner1.yRange  = [0, 0]
forceCorner1.zRange  = [0, 0]
forceCorner1.value   = -0.0005

# Last boundary condition sets the (1,1,0) corner to be half the line force
forceCorner2.type    = HypOptLib.BoundaryConditionType.LOAD
forceCorner2.degreesOfFreedom = {2}
forceCorner2.xRange  = [2, 2]
forceCorner2.yRange  = [1, 1]
forceCorner2.zRange  = [0, 0]
forceCorner2.value   = -0.0005

# Apply boundary conditions
solver.setBoundaryConditions([fixedPoints, forceCentre, forceCorner1, forceCorner2])

######################################################################
# Setup Hyperoptimization parameters.
solver.setSavePath("mbb_beam.h5")
solver.setTargetTemperature(0.01)
solver.setTimestep(0.01)
solver.setMaximumIterations(6000)

######################################################################
# Run the simulation
saveRange = [0,6000]
solver.newRun(saveRange)
