
=======================================
Python Wrapper
=======================================

This section describes the functionality exposed to Python in HypOptLib.

Note that all parameter types here are shown as C++ variables. The Python wrapper auto-translates
these from Python types. For example, a Python *list of ints* will be translated into a C++ 
*std::vector<uint32_t>** automatically.

Contents
======================================

* :ref:`Structs and Enums <Structs and Enums>`
    * :ref:`BoundaryConditionType <BoundaryConditionType>`
    * :ref:`verbosity <verbosity>`
    * :ref:`BoundaryCondition <BoundaryCondition>`
    * :ref:`DomainCoordinates <DomainCoordinates>`

* :ref:`HypOptLib <HypOptLib_header>`

    .. hlist::
        :columns: 2

        * :ref:`Basic Functionality <Basic Functionality>`
            * :ref:`HypOptLib <HypOptLib_function>`
            * :ref:`newRun <newRun>`
            * :ref:`restartRun <restartRun>`

        * :ref:`Advanced Functionality <Advanced Functionality>`
            * :ref:`generateRandomInitialConditionsFile <generateRandomInitialConditionsFile>`

        * :ref:`Required Settings <Required Settings>`
            * :ref:`setGridProperties <setGridProperties>`
            * :ref:`setBoundaryConditions <setBoundaryConditions>`

        * :ref:`Recommneded Settings <Recommneded Settings>`
            * :ref:`setSavePath <setSavePath>`
            * :ref:`setTargetTemperature <setTargetTemperature>`
            * :ref:`setTimestep <setTimestep>`
            * :ref:`setNoseHooverChainOrder <setNoseHooverChainOrder>`
            * :ref:`setMaximumIterations <setMaximumIterations>`

        * :ref:`Optional Settings <Optional Settings>`
            * :ref:`setPenalty <setPenalty>`
            * :ref:`setMinimumFilterRadius <setMinimumFilterRadius>`
            * :ref:`setVolumeFraction <setVolumeFraction>`
            * :ref:`setRandomStartingValues <setRandomStartingValues>`
            * :ref:`setSaveHamiltonian <setSaveHamiltonian>`
            * :ref:`setMaxSimulationTime <setMaxSimulationTime>`
            * :ref:`enableVariableTimestep <enableVariableTimestep>`
            * :ref:`loadInitialConditionsFromFile <loadInitialConditionsFromFile>`
            * :ref:`restartDoesntSupportCustomMesh <restartDoesntSupportCustomMesh>`
            * :ref:`setMaximumFeaSolverIterations <setMaximumFeaSolverIterations>`

.. _Structs and Enums:

Structs and Enums
=======================================

These objects are used to define the physical characteristics of the problem: the mesh and the boundary
conditions.

.. _BoundaryConditionType:

.. doxygenenum:: BoundaryConditionType
   :project: HypOptLib

.. _verbosity:

.. doxygenenum:: verbosity
   :project: HypOptLib

.. _BoundaryCondition:

.. doxygenstruct:: BoundaryCondition
   :project: HypOptLib
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. _DomainCoordinates:

.. doxygenstruct:: DomainCoordinates
   :project: HypOptLib
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. _HypOptLib_header:

HypOptLib
=======================================

.. doxygenclass:: HypOptLib
   :project: HypOptLib

.. _Basic Functionality:

Basic Functionality
---------------------------------------

The following functions are used to initialize and start the hyperoptimization algorithm.


.. _HypOptLib_function:
.. doxygenfunction:: HypOptLib
    :project: HypOptLib

.. _newRun:

.. doxygenfunction:: newRun
    :project: HypOptLib

.. _restartRun:

.. doxygenfunction:: restartRun
    :project: HypOptLib


.. _Advanced Functionality:

Advanced Functionality
---------------------------------------

For more advanced hyperoptimization usage or debugging.

.. _generateRandomInitialConditionsFile:

.. doxygenfunction:: generateRandomInitialConditionsFile
    :project: HypOptLib


.. _Required Settings:

Required Settings
---------------------------------------

These functions are required in order to define the topology optimization problem.

.. _setGridProperties:

.. doxygenfunction:: setGridProperties
    :project: HypOptLib

.. _setBoundaryConditions:

.. doxygenfunction:: setBoundaryConditions
    :project: HypOptLib

.. _Recommneded Settings:

Recommneded Settings
---------------------------------------

These functions are recommended for basic usage. Each will use default values if not set,
but are commonly set to custom values.

.. _setSavePath:

.. doxygenfunction:: setSavePath
    :project: HypOptLib

.. _setTargetTemperature:

.. doxygenfunction:: setTargetTemperature
    :project: HypOptLib

.. _setTimestep:

.. doxygenfunction:: setTimestep
    :project: HypOptLib

.. _setNoseHooverChainOrder:

.. doxygenfunction:: setNoseHooverChainOrder
    :project: HypOptLib

.. _setMaximumIterations:

.. doxygenfunction:: setMaximumIterations
    :project: HypOptLib


.. _Optional Settings:

Optional Settings
---------------------------------------

These functions help specify specific parameters for hyperoptimization runs. If not specified,
each parameter will be set to a default value.

.. _setPenalty:

.. doxygenfunction:: setPenalty
    :project: HypOptLib

.. _setMinimumFilterRadius:

.. doxygenfunction:: setMinimumFilterRadius
    :project: HypOptLib

.. _setVolumeFraction:

.. doxygenfunction:: setVolumeFraction
    :project: HypOptLib

.. _setRandomStartingValues:

.. doxygenfunction:: setRandomStartingValues
    :project: HypOptLib

.. _setSaveHamiltonian:

.. doxygenfunction:: setSaveHamiltonian
    :project: HypOptLib

.. _setMaxSimulationTime:

.. doxygenfunction:: setMaxSimulationTime
    :project: HypOptLib

.. _enableVariableTimestep:

.. doxygenfunction:: HypOptLib::enableVariableTimestep
    :project: HypOptLib

.. _loadInitialConditionsFromFile:

.. doxygenfunction:: loadInitialConditionsFromFile
    :project: HypOptLib

.. _restartDoesntSupportCustomMesh:

.. doxygenfunction:: restartDoesntSupportCustomMesh
    :project: HypOptLib

.. _setMaximumFeaSolverIterations:

.. doxygenfunction:: setMaximumFeaSolverIterations
    :project: HypOptLib
