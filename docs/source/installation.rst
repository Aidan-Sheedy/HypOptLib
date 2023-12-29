========================
Installation and Setup
========================

HypOptLib is designed to run in Linux environments. It has been tested mainly
on WSL and Ubuntu, but should in theory support any Linux distribution with the
appropriate requirements.

HypOptLib can be setup in two main ways. A conveniance script `macros.sh` is
provided which can automatically install any prerequisites along with other
functionality. Alternatively, the prerequisites and their recommended method of
installation are also provided here. For full documentation of the macros script
see this link 

.. todo::

    include link.

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
of the following requirements:

+------------+--------------------------------------------------+
| Dependency | Command (Ubuntu)                                 |
+============+==================================================+
| cmake      | sudo apt install cmake                           |
+------------+--------------------------------------------------+
| make       | sudo apt install make                            |
+------------+--------------------------------------------------+
| mpi        | sudo apt install mpich                           |
+------------+--------------------------------------------------+
| Pybind11   | pip3 install "pybind11[global]"                  |
+------------+--------------------------------------------------+
| HDF5       | sudo apt install libhdf5-serial-dev              |
+------------+--------------------------------------------------+
| BLAS       | sudo apt-get install libblas-dev liblapack-dev   |
+------------+--------------------------------------------------+

Finally install PETSc, following the instructions here: https://petsc.org/release/install/download/.

When configuring PETSc, you must add the `--download-hdf5` tag:

.. code-block:: bash

    ./configure --download-hdf5

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
A few examples are provided in `run/main.py`, but a basic script works as follows:

.. code-block:: python

    #!/usr/bin/env python3

    import HypOptLib
    solver = HypOptLib.HypOptLib()

    # Set up solver settings
    solver.setTargetTemperature(0.1)
    solver.setTimestep(0.001)
    solver.setMaximumIterations(1000)

    saveRange = [900, 1000]
    gridDeimensions = [32, 16, 16]

    # Start Simulation
    solver.newRun( saveRange, gridDeimensions )

This basic script can then ammended with all the specific settings applicable to
the desired simulation. Full documentation is provided here. 

.. todo::

    Link this to the full documentation.






