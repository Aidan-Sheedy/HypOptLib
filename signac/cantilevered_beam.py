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

class CantileveredBeam:

    solver = HypOptLib.HypOptLib()
    domain = HypOptLib.DomainCoordinates()

    ######################################################################
    # Setup Hyperoptimization parameters.

    def __init_boundary_conditions(self):
        fixedPoints = HypOptLib.BoundaryCondition()
        forceCentre = HypOptLib.BoundaryCondition()
        forceCorner1  = HypOptLib.BoundaryCondition()
        forceCorner2  = HypOptLib.BoundaryCondition()

        ######################################################################
        # Setup domain. This describes a rectangular prism of dimensions 2x1x1, with cubic voxels.
        self.domain.xMinimum = 0
        self.domain.xMaximum = 2
        self.domain.yMinimum = 0
        self.domain.yMaximum = 1
        self.domain.zMinimum = 0
        self.domain.zMaximum = 1

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

        self.solver.setBoundaryConditions([fixedPoints, forceCentre, forceCorner1, forceCorner2])

        return


    def __init__(self, grid, timestep, maxIterations, temperature, saveFreq, filePath, volfrac=0.5, saveHamiltonian=True):

        initial_conditions_path = filePath[:-3] + "_initial_conditions" + filePath[-3:]
        self.__init_boundary_conditions()
        self.solver.setTargetTemperature(temperature)
        self.solver.setGridProperties(grid, self.domain)
        self.solver.setTimestep(timestep)
        self.solver.setMaximumIterations(maxIterations)
        self.solver.setVolumeFraction(volfrac)
        self.solver.setSaveHamiltonian(saveHamiltonian)
        self.solver.setSaveFrequency(saveFreq)

        self.solver.generateRandomInitialConditionsFile(grid, initial_conditions_path)
        self.solver.loadInitialConditionsFromFile(initial_conditions_path)

        self.solver.setSavePath(filePath)

    def startRun(self, saveRange):
        self.solver.newRun(saveRange)

