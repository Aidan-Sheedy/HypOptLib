=======================================
Library Documentation
=======================================



Class Structure
=======================================

HypOptLib is built around the original TopOpt code, replacing the MMA algorithm with a new Hyperoptimization class.
The Hyperoptimization class uses all Petsc types and functions, and has full parallelization support. The call
structure is shown in :numref:`Fig1`, with the Pyhon library itself shown in blue. The library takes in parameters
defined in Python, and instantiates all the classes required to run the design loop. The loop itself is then written
and run in the Hyperoptimization class itself. Wrappers are provided for the Lagrange multipliers, filtering, and
sensitivity calculations, shown in green. These define the functions specific to whatever problem is being solved. In
the package provided, they define a compliance-optimized cantilivered beam. Some of this functionality is provided by
the original TopOpt classes, shown highlighted in red in :numref:`Fig1`.

There are three supporting classes unrelated to the solver itself. TopOpt is mostly used to initialize the finite
element analysis solver and vectors, while the original use for storing all parameters is no longer used. The
FileManager class provides all HDF5 file operations, and reads or writes the final outputs as required. Finally, a
utility class PetscExtensions is provided to manage some parallelization functions not present in Petsc itself.

.. _Fig1:

.. figure:: ../figures/software_diagram.svg
    :width: 800

    Diagram indicating the class structure in HypOptLib. Classes at the top call classes below them.
    Smaller classes on the right indicate supporting classes. HypOptLib is the library exposed via Pybind11, and
    instantiates all classes. Hyperoptimization then calls the green wrapper classes, which provide lower level
    functionality.

How to Modify HypOptLib
=======================================

If you want to update HypOptLib to apply hyperoptimization to a new problem, there are two ways depending on the
problem.

Topology optimization
---------------------------------------

For topology optimization, some properties are exposed to the Python wrapper, such as:

 * Volume Fraction,
 * Filter radius,
 * Penalty,
 * Grid dimensions.

However, the rest of the problem definition is still written in C++. For an arbitrary new Topology optimization
problem, the following classes may need to be modified:

 1. **LinearElasticity** must be modified for all boundary conditions.
 2. **LagrangeMultiplier** must be used as a base class for any new overwritten Lagrange Multiplier calculations.
 3. **FilterWrapper** must be used as a base class for any new filters. The filter *must* be linear. **Hazhir is this true?**
 4. **SensitivityWrapper** must be used as a base class for sensitivity and objective function calculations. 

The way the last three classes are written, any new implementation can create a new derived class which implements all
the virtual functions provided in the base wrapper classes. This allows for simpler code that does not need to go and
re-write old functions. For example, if all the original classes can be reused for a given problem, but the
LagrangeMultiplier needs to be redone, then a new **ExampleLagrangeMultiplier** class can be created which uses the
original **LagrangeMultiplier** class as a base, but overwrites the LagrangeMultiplier::computeLagrangeMultiplier
function. This can then be passed to the Hyperoptimization when initializing, and the new
ExampleLagrangeMultiplier::computeLagrangeMultiplier function will be used without having to change the original
Hyperoptimization code.

Generic optimization problems
---------------------------------------

For any other optimization problem, the above section still holds true. However, it is possible (if not likely) that
the HypOptLib class will need to be rewritten. For ease of use, it is recommended that as many parameters be exposed to
the Python wrapper if possible if changes are being made to the HypOptLib class, to avoid having to re-compile the
library whenver changes are needed.


Additional Information
=======================================

There are two main additional features included with HypOptLib: the macros script and analysis scripts. These are both
provided as optional tools to make running and using HypOptLib easier, but are not required in the slightest.

Analysis Scripts
---------------------------------------

Several Python scripts are provided for basic data analysis. They are as follows:

.. warning::
    **TODO** Fill out the rest of this section

macros.sh
---------------------------------------

The macros script is intended to take all the command line operations usually necessary for building, compiling,
installing, etc. and reduce them to one script. For example, with cmake a fresh build will require making a build
folder, running cmake, and finally running make:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make

While this might not seem difficult, when actively working on a project this can get tiresome - especially when adding
lots of files, or debugging makefile or cmake issues. As such, actions such as this have been reduced to a single
command line option in the build.sh script.

The full capabilities of the script are as follows:

.. warning::
    **TODO** Layout the options for the macros script
