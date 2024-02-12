========================
Installation and Setup
========================

HypOptLib is designed to run in Linux environments. It has been tested mainly
on WSL and Ubuntu, but should in theory support any Linux distribution with the
appropriate requirements.

HypOptLib can be setup in two main ways. A conveniance script `macros.sh` is
provided which can automatically install any prerequisites along with other
functionality. Alternatively, the prerequisites and their recommended method of
installation are also provided here. The full documentation of suported macros is
provided :doc:`here </apidocs/macros>`.

Automatic Install
========================

1. Clone the HypOptLib repository from https://github.com/Aidan-Sheedy/Hyperoptimization_using_Petsc.git.

2. Run the automatic install macro:

    .. code-block:: bash

        ./macros.sh setup

    .. note:: 
        
        You may need to change the executable permissions first:
        
        .. code-block:: bash

            chmod +rwx ./macros.sh

3. Restart the console

Manual Install
========================

For manual install, clone the HypOptLib repository as above, then install each
of the following requirements. First, you may have to run

.. code-block:: bash

    sudo apt update

+-------------+--------------------------------------------------+
| Dependency  | Command (Ubuntu)                                 |
+=============-+==================================================+
| cmake       | sudo apt install cmake                           |
+-------------+--------------------------------------------------+
| make        | sudo apt install make                            |
+-------------+--------------------------------------------------+
| mpi         | sudo apt install mpich                           |
+-------------+--------------------------------------------------+
| python pip3 | sudo apt install python3-pip                     |
+-------------+--------------------------------------------------+
| Pybind11    | pip3 install "pybind11[global]"                  |
+-------------+--------------------------------------------------+
| HDF5        | sudo apt install libhdf5-serial-dev              |
+-------------+--------------------------------------------------+
| BLAS        | sudo apt-get install libblas-dev liblapack-dev   |
+-------------+--------------------------------------------------+

Finally install PETSc, following the instructions here: https://petsc.org/release/install/download/.

When configuring PETSc, you must add the `--download-hdf5` tag:

.. code-block:: bash

    ./configure --download-hdf5

Some operating systems or environments may also need additional configure flags. If the compilation
fails, check the `PETSc install reference <https://petsc.org/release/install/install/>` for help.

.. note::

    Some package managers carry PETSc, but make sure that the version is 3.17 or newer.


========================
Building and Running
========================

They macros script again provides automatic building, but this can still be done
manually.

Manual Buildings
========================

1. Make a folder called `build` in the HypOptLib directory.

2. In that folder, run:
    .. code-block:: bash

        cmake ..

3. Then, run:
    .. code-block:: bash

        make

4. The output will be built in the `run` folder alongside `main.py`. Alternatively, 


Automatic Building
========================

To build automatically, use the following macro.sh command:

.. code-block:: bash

    ./macros.sh build [clean/all]

By default, the `all` build option will be used, which builds cmake and make commands. If
the cmake output is already complete, only make will be run. The `clean` build option will
clear all cmake and make outputs and objects.

Running HypOptLib
========================

The HypOptLib library can be imported just like any other Python library. The library binary
can either be added to the PATH, or can simply be in the same directory as the Python script.

To run a python script with MPI, use the command:

.. code-block:: bash

    mpiexec -n x Python3 ./pythonScript.py

where x is the number of cores to run on. It is recommended to use an even number of cores.

A few example scripts are provided in `examples/`, as well as a barebones script in `run/main.py`,
but a basic script works as follows:

.. code-block:: python

    #!/usr/bin/env python3

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

    gridDeimensions = [32, 16, 16]
    solver.setGridProperties(gridDeimensions, domain)

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

    solver.setBoundaryConditions([fixedPoints, forceCentre, forceCorner1, forceCorner2])

    # Set up solver settings
    solver.setTargetTemperature(0.1)
    solver.setTimestep(0.001)
    solver.setMaximumIterations(1000)
    saveRange = [900, 1000]

    # Start Simulation
    solver.newRun(saveRange)

This basic script can then ammended with all the specific settings applicable to
the desired simulation. Full documentation is provided :doc:`here </apidocs/pybind11>`.
