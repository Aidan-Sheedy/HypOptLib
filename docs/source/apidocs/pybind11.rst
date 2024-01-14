
=======================================
HypOptLib - Python Wrapper
=======================================

The Python wrapper exposes the HypOptLib class to Python. All public functions are exposed,
as described below.

.. note::

    All parameter types here are shown as C++ variables. The Python wrapper auto-translates these
    from Python types. For example, a Python *list of ints* will be translated into a C++ 
    *std::vector<uint32_t>** automatically.

Basic Functionality
=======================================

The following functions are used to initialize and start the hyperoptimization algorithm.

.. doxygenfunction:: HypOptLib
    :project: HypOptLib

.. doxygenfunction:: newRun
    :project: HypOptLib

.. doxygenfunction:: restartRun
    :project: HypOptLib

Advanced Functionality
=======================================

For more advanced hyperoptimization usage or debugging.

.. doxygenfunction:: generateRandomInitialConditionsFile
    :project: HypOptLib

Optional Settings
=======================================

These functions help specify specific parameters for hyperoptimization runs. If not specified,
each parameter will be set to a default value.

.. doxygenfunction:: setSavePath
    :project: HypOptLib

.. doxygenfunction:: setTargetTemperature
    :project: HypOptLib

.. doxygenfunction:: setTimestep
    :project: HypOptLib

.. doxygenfunction:: setNoseHooverChainOrder
    :project: HypOptLib

.. doxygenfunction:: setMaximumIterations
    :project: HypOptLib

.. doxygenfunction:: setPenalty
    :project: HypOptLib

.. doxygenfunction:: setMinimumFilterRadius
    :project: HypOptLib

.. doxygenfunction:: setVolumeFraction
    :project: HypOptLib

.. doxygenfunction:: setRandomStartingValues
    :project: HypOptLib

.. doxygenfunction:: setSaveHamiltonian
    :project: HypOptLib

.. doxygenfunction:: setMaxSimulationTime
    :project: HypOptLib

.. doxygenfunction:: HypOptLib::enableVariableTimestep
    :project: HypOptLib

.. doxygenfunction:: loadInitialConditionsFromFile
    :project: HypOptLib



.. .. admonition:: HypOptLib()

..     Empty constructor, no initialization done here.

.. .. admonition:: newRun(iterationSaveRange: "list of ints", gridDimensions: "list of ints")

..     Starts a fresh run with the provided parameters.

..     **Parameters:**

..         * **iterationSaveRange** Range of iterations to save. Does not need to include the final
..         iteration to support restarting, this is saved regardless. Can be set to [0,0] to disable
..         saving any iterations.

..         * **gridDimensions** Dimensions of the cells in the grid. The final mesh will have be one
..         higher in each dimension. The dimenions should each be divisible by 2 three times.

..     **Throws:** HypOptException.

..     **Returns:** 0 on success, or error. Only used to align with Petsc exception handling.

.. .. admonition:: restartRun(restartPath: str, iterationSaveRange: "list of ints")

..     Restarts a simualtion from the provided file path. All options are parsed from the metadata in
..     the restart file. This means that any optional settings provided will be ignored.

..     **Parameters:**

..         * **restartPath** Path to the hdf5 file to restart.

..         * **iterationSaveRange** Range of iterations to save. Does not need to include the final
..         iteration to support restarting, this is saved regardless. Can be set to [0,0] to disable
..         saving any iterations. Iterations restart at 1, not from the last iteration of the provided
..         file.

..     **Throws:** HypOptException.

..     **Returns:** 0 on success, or error. Only used to align with Petsc exception handling.

.. Advanced Functionality
.. =======================================

.. .. admonition:: generateRandomInitialConditionsFile(gridDimensions: "list of ints", filePath: str)

..     Generates a file with randomized initial velocities and positions. This file can then be passed
..     as an optional parameter to ensure the same initial conditions are used accross multiple runs.

..     **Parameters:**

..         * **gridDimensions** Dimensions of the cells in the grid. The final mesh will have be one
..         higher in each dimension. The dimenions should each be divisible by 2 three times.

..         * **filePath** Name and path to save the output to.

..     **Throws:** HypOptException.

