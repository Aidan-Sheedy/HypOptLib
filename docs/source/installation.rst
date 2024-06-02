++++++++++++++++++++++++
Getting Started
++++++++++++++++++++++++

These sections describe downloading, building, and running HypOptLib.
This should act as a basic quick-start tutorial.

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

This installation method can be useful for personal systems, but may not be ideal as it
won't account for any idiosyncrasies. However, it can also be useful to look through the
bash script to see which dependencies are being installed.

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
of the following requirements. You may have to run

.. code-block:: bash

    sudo apt update

before installing dependencies.

+-------------+--------------------------------------------------+
| Dependency  | Command (Ubuntu)                                 |
+=============+==================================================+
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
fails, check the `PETSc install reference <https://petsc.org/release/install/install/>`_ for help.

.. note::

    Some package managers carry PETSc, but make sure that the version is 3.17 or newer.


========================
Building and Running
========================

The macros script again provides automatic building, but this can still be done
manually.

In general, a few items may need to be set manually to ensure the build works:

1. In **HypOptLib/CMakeLists.txt**, the two lines setting **PETSC_DIR** and **PETSC_ARCH**
   may need to be uncommented and set to point to your local PETSc install location.

2. In the same file, the two lines setting the **PythonInterp** and **PythonLibs** versions
   may be needed to help CMake find the correct Python version.

3. At the end of the file, you may need to uncomment the line 

   .. code-block:: bash

        # list( APPEND CMAKE_INSTALL_RPATH ${PETSC_DIR}/lib )

   if petsc is not being found or if the library is not being imported to python after building
   with make install. Only do this if you actually want to add the PETSc library to the executable
   rpath.

Manual Buildings
========================

1. Make a folder called `build` in the HypOptLib directory.

2. In that folder, run:

    .. code-block:: bash

        cmake ..

3. Then, run:

    .. code-block:: bash

        make

4. The output binary will be built in the `run` folder alongside `main.py`. You can import HypOptLib
   in any Python file that has this binary in the same runtime directory.

Make Install
------------------------

To install HypOptLib system or user-wide, you can run :code:`make install` to install the library
to the location determined by the cmake command.

The default install location is to **${CMAKE_INSTALL_PREFIX}/lib**, where for my WSL setup this resolves
to **/usr/local/lib**. This is probably not the ideal location for an install, so there are a few cmake
command line options to help set a better install location:

 * -DCUSTOM_INSTALL_PREFIX=<prefix> will set the install path to **${prefix}/lib** instead of
   **${CMAKE_INSTALL_PREFIX}/lib**

 * -DFORCE_INSTALL_PATH=<path> will set the install path to **${CMAKE_INSTALL_PREFIX}/<path>** instead of
   **${CMAKE_INSTALL_PREFIX}/lib**

 * -DInstallPythonSysPath=ON will use the Python_SITELIB as the libray install directory. This is **not**
   recommended in general, especially on clusters or shared systems. However it is a quick and dirty way to
   install the library globally on personal systems.

Of course, it these options may not be enough for all circumstances, so it may still be necessary to modify
the cmake files for your own specific environment.

Automatic Building
========================

To build automatically, use the following macro.sh command:

.. code-block:: bash

    ./macros.sh build [clean/all]

By default, the `all` build option will be used, which builds cmake and make commands. If
the cmake output is already complete, only make will be run. The `clean` build option will
clear all cmake and make outputs and objects. This should be equivalent to following the manual
installation without running :code:`make install`.

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

========================
Building Documentation
========================

If for some reason you need to build the documentation locally, follow these steps.

First, install these dependencies:

+-------------+--------------------------------------------------+
| Dependency  | Command (Ubuntu)                                 |
+=============+==================================================+
| doxygen     | sudo apt install doxygen                         |
+-------------+--------------------------------------------------+
| sphinx      | sudo apt install sphinx                          |
+-------------+--------------------------------------------------+
| breathe     | pip3 install breathe                             |
+-------------+--------------------------------------------------+
| rtd_theme   | pip3 install sphinx_rtd_theme                    |
+-------------+--------------------------------------------------+

Then, run

.. code-block:: bash

    ./macros.sh build_docs

to compile the documentation. The output will be available in **docs/build/html/index.html**.
