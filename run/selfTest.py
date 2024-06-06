#!/usr/bin/env python3

'''
This file implements self-testing functionality for HypOptLib. It is
intended to confirm if an installation has been successful.

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

import numpy as np
import contextlib
import HypOptLib
import h5py
import sys
import os

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


    def __init__(self, grid, timestep, maxIterations, temperature, saveFreq, filePath, initialConditions="", volfrac=0.5, saveHamiltonian=False):

        self.__init_boundary_conditions()
        self.solver.setTargetTemperature(temperature)
        self.solver.setGridProperties(grid, self.domain)
        self.solver.setTimestep(timestep)
        self.solver.setMaximumIterations(maxIterations)
        self.solver.setVolumeFraction(volfrac)
        self.solver.setSaveHamiltonian(saveHamiltonian)
        self.solver.setSaveFrequency(saveFreq)

        if ("random" == initialConditions):
            initial_conditions_path = filePath[:-3] + "_initial_conditions" + filePath[-3:]
            self.solver.generateRandomInitialConditionsFile(grid, initial_conditions_path)
            self.solver.loadInitialConditionsFromFile(initial_conditions_path)
        elif ("" != initialConditions):
            self.solver.loadInitialConditionsFromFile(initialConditions)
        else:
            self.solver.setRandomStartingValues(False)

        self.solver.setSavePath(filePath)

    def startRun(self, saveRange):
        self.solver.newRun(saveRange)


class Tester:

    def basicSanityTest(self):
        '''
        Tests a simple uniform case with 0 temperature. 
        '''
        baseline="../selfTestFiles/basicSantiy_t0_16x8x8_200itr.h5"

        dimensions      = [16, 8, 8]
        timestep        = 0.1
        maxIterations   = 200
        temperature     = 0
        saveFile        = "basicSantiyTestResults.h5"

        beam = CantileveredBeam(dimensions,
                                timestep,
                                maxIterations,
                                temperature,
                                maxIterations,
                                saveFile)

        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                beam.startRun([0,0])

        return self.compareResults(saveFile, baseline)

    def compareResults(self, test, baseline, threshold=0.00001, compareCompliance=False):
        results = [threshold, True, True, True]

        testFile = h5py.File(test)
        baseFile = h5py.File(baseline)
        testFinalPosition   = np.asarray(testFile["/Dataset/Final Position"])
        baseFinalPosition   = np.asarray(baseFile["/Dataset/Final Position"])
        testTemperature     = np.asarray(testFile["/Dataset/Temperature"])
        baseTemperature     = np.asarray(baseFile["/Dataset/Temperature"])
        if (compareCompliance):
            testCompliance  = np.asarray(testFile["/Dataset/Compliance"])
            baseCompliance  = np.asarray(baseFile["/Dataset/Compliance"])
        testFile.close()
        baseFile.close()

        if (not np.allclose(testFinalPosition, baseFinalPosition, rtol=threshold, atol=0)):
            results[1] = False

        if (not np.allclose(testTemperature, baseTemperature, rtol=threshold, atol=0)):
            results[2] = False

        if (compareCompliance):
            if (not np.allclose(testCompliance, baseCompliance, rtol=threshold, atol=0)):
                results[3] = False

        return results

    def parseResults(self, resultsArray, compareCompliance=False):
        percentage = resultsArray[0]*100
        if (resultsArray[1]):
            print("Positions:\tPASS - Match within " + str(percentage) + "%")
        else:
            print("Positions:\tFAIL - Not within " + str(percentage) + "%")

        if (resultsArray[2]):
            print("Temperatures:\tPASS - Match within " + str(percentage) + "%")
        else:
            print("Temperatures:\tFAIL - Not within " + str(percentage) + "%")
 
        if (compareCompliance):
            if (resultsArray[3]):
                print("Compliances:\tPASS - Match within " + str(percentage) + "%")
            else:
                print("Compliances:\tFAIL - Not within " + str(percentage) + "%")

        return not (False in resultsArray)


def main():
    status = True
    tester = Tester()

    print("\n============  Start Self Testing  ============")
    print("------------  Test: Basic Sanity  ============")

    basicSanityResults = tester.basicSanityTest()

    print("============ Self Testing Results ============")
    print("------------  Test: Basic Sanity  ============")
    status = tester.parseResults(basicSanityResults)

    print("\n++++++++++++++++++++ NOTE ++++++++++++++++++++")
    print("This self-test is currently limited to running one basic")
    print("instance and checking against expected results. The following")
    print("tests are not checked, and may or may not be added in the future.")
    print(" - Parallel MPI tests")
    print(" - Randomized starting conditions")
    print(" - Negative testing (ensuring changing variables changes results)")
    print(" - Qualitative testing (Ensuring that simulation behaviour is as expected)")
    print(" - Comparison to analytical results (this may not be feasible)")
    print("\nNo promises on if or when these will be implemented :)")
    print("\nAlso, the current test compares against the values generated")
    print("on the author's personal computer from the commit hash")
    print("18d79b40bdea0da6ceaed8a0b63d2a1110aced22. The author is aware")
    print("of the shortcomings of this approach.")

    return 0 if True == status else -1

if __name__ == "__main__":
    sys.exit(main())
