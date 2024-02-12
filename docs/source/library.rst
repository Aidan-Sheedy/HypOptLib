=======================================
Library Documentation
=======================================


Class Structure
=======================================

HypOptLib is built around the original TopOpt code, replacing the MMA algorithm with a new Hyperoptimization class.
The Hyperoptimization class uses all Petsc types and functions, and has full parallelization support. The call
structure is shown in :numref:`Fig1`, and follows these steps:

    1. The problem is defined in a Python script. This includes boundary conditions and the mesh, the temperature,
       timestep, iteration count, and a variety of other optional settings.

    2. The Python file is interpreted by a Pybind11 wrapper, which wraps the HypOptLib class. This class contains basic
       settings functions, such as options for meshes and boundary conditions, temperature, timestep, number of iterations, etc.

    3. HypOptLib instantiates the hyperoptimization problem based on the settings provided by the problem file, and when indicated
       starts the hyperoptimization simulation.

    4. Each iteration, Hyperoptimization calls a number of wrappers for Lagrange Multiplier, Filter, and Sensitivity calculations.
       These then either call implementation classes (provided by TopOpt by default), or implement the calculation directly.

    5. At the end of each iteration, and at the end of the simulation, data is saved by the FileManager class. Files are saved in the
       HDF5 format.

    6. Finally, all classes from Hyperoptimization and down utilize calls to PETSc. Some small utility functions are also provided
       by an extension class.

It should be noted that three classes are used from the original TopOpt library: TopOpt, Filter, and Linear Elasticity.

.. _Fig1:

.. figure:: ../figures/software_diagram.svg
    :width: 800

    Diagram indicating the class structure in HypOptLib. Classes at the top call classes below them.
    HypOptLib is the library exposed via Pybind11, and instantiates all classes. Hyperoptimization then
    calls the green wrapper classes, which provide lower level functionality. All classes below and including
    Hyperoptimization make calls to the PETSc library for computations and solvers.

How to Modify HypOptLib
=======================================

If you want to update HypOptLib to apply hyperoptimization to a new problem, there are two ways depending on the
problem.

Topology optimization
---------------------------------------

For topology optimization, some properties are exposed to the Python wrapper, such as:

 * Volume fraction,
 * Filter radius,
 * Penalty,
 * Grid dimensions,
 * Boundary conditions,
 * Temperature,
 * Nose Hoover chain order.

This should cover most generic topology optimization problems, but a few small features are not yet implemented:

 * Void points or arbitrary-density fixed points,
 * Non-rectilinear boundary conditions,
    * This can be done manually in Python by generating arrays of single-point boundary conditions or something similar.
      However, this might not be ideal and users may want to expand on the exisitng boundary condition options.
 * Advanced initial condition files,
    * This can be done by generating an HDF5 initial conditions file that matches the data format used by the library.
      However, there is no support for standard 3D modelling data formats as initial condition inputs, such as stl.
 * non-rectilinear meshes.

Additionally, changes to the filtering, Lagrange Multipliers, or Linear Elasticity calculations must be done in C++,
as the have not been exposed to Python.

Generic optimization problems
---------------------------------------

For any other optimization problem, the above section still holds true. However, new sensitivity, filter, and Lagrange Multiplier
calcualtions will probably be necessary. However, there is no need to overwrite the existing classes. Instead, derived classes of
each type can be created which inherit the original wrapper functions. So long as the derived classes have new implementations of
the virtual functions provided in the wrapper, the Hyperoptimization class will call the derived implementations instead of the
wrapper.

For example, if all the original classes can be reused for a given problem but the LagrangeMultiplier needs to be redone, then a
new **ExampleLagrangeMultiplier** class can be created. This will inherit the original **LagrangeMultiplier** class, but will
overwrite **LagrangeMultiplier::computeLagrangeMultiplier**. **ExampleLagrangeMultiplier** can then be passed to **Hyperoptimization**
when initializing, and the new **ExampleLagrangeMultiplier::computeLagrangeMultiplier** function will be used without having to change
the original Hyperoptimization code.

The following files can be changed based on the requirements for the new optimization problem:

 1. **LinearElasticity** must be modified for all boundary conditions.
 2. **LagrangeMultiplier** must be used as a base class for any new overwritten Lagrange Multiplier calculations.
 3. **FilterWrapper** must be used as a base class for any new filters. The filter *must* be linear. **Hazhir is this true?**
 4. **SensitivityWrapper** must be used as a base class for sensitivity and objective function calculations.

Finally, if for any reason HypOptLib needs to be edited (for example to add additional parameters to save), it is recommended that
as much is exposed to the Python wrapper as possible to avoid re-compiling the library.
