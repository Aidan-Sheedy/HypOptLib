#!/usr/bin/env python3

'''
This file shows a stripped-down hyperoptimization example as a starting point.

See the examples folder for specific examples, and the documentation for full
api information.

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
#
# First boundary condition fixes the x=0 plane
fixedPoints.type    = HypOptLib.BoundaryConditionType.FIXED_POINT
fixedPoints.xRange  = [0, 0]
fixedPoints.yRange  = [0, 1]
fixedPoints.zRange  = [0, 1]
fixedPoints.degreesOfFreedom = {0, 1, 2}
fixedPoints.value   = 0

# Second boundary condition sets a line force at X=1, Z=0.5, in the Z DOF
forceCentre.type    = HypOptLib.BoundaryConditionType.LOAD
forceCentre.xRange  = [2, 2]
forceCentre.yRange  = [0, 1]
forceCentre.zRange  = [0.5, 0.5]
forceCentre.degreesOfFreedom = {2}
forceCentre.value   = -0.001

# Third boundary condition sets the (1,0,0.5) corner to be half the line force
forceCorner1.type    = HypOptLib.BoundaryConditionType.LOAD
forceCorner1.degreesOfFreedom = {2}
forceCorner1.xRange  = [2, 2]
forceCorner1.yRange  = [0, 0]
forceCorner1.zRange  = [0.5, 0.5]
forceCorner1.value   = -0.0005

# Last boundary condition sets the (1,1,0.5) corner to be half the line force
forceCorner2.type    = HypOptLib.BoundaryConditionType.LOAD
forceCorner2.degreesOfFreedom = {2}
forceCorner2.xRange  = [2, 2]
forceCorner2.yRange  = [1, 1]
forceCorner2.zRange  = [0.5, 0.5]
forceCorner2.value   = -0.0005

solver.setBoundaryConditions([fixedPoints, forceCentre, forceCorner1, forceCorner2])

######################################################################
# Setup Hyperoptimization parameters.

solver.setTimestep(0.01)
solver.setMaximumIterations(100)
saveRange = [0, 100]
solver.newRun(saveRange)
