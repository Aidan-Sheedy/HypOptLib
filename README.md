# HypOptLib

HypOptLib is a Petsc based semi-parallel implementation of hyperoptimization with a Python wrapper. The package is a
fork of the topology optimization framework
[TopOpt](https://www.topopt.mek.dtu.dk/apps-and-software/large-scale-topology-optimization-code-using-petsc) by Aage
et. al., along with some borrowed elements from the derivative [TopOptLib](https://doi.org/10.1007/s00158-021-03018-7) by Smit et. al.

You can find the full documentation at https://aidan-sheedy.github.io/HypOptLib/.

## Credit

**Code** Aidan Sheedy

**Hyperoptimization Algorithm** Hazhir Aliahmadi

**Advisor** Greg van Anders

## Borrowed Material

**TopOpt** [Niels Aage, Erik Andreassen, Boyan Stefanov Lazarov](https://www.topopt.mek.dtu.dk/apps-and-software/large-scale-topology-optimization-code-using-petsc)

**TopOptLib** [Thijs Smit, Niels Aage, Stephen J. Ferguson, Benedikt Halgason](https://doi.org/10.1007/s00158-021-03018-7)

## Is HypOptLib right for you?

HypOptLib is intended to demonstrate the possibilities of hyperoptimization, to show off what it can do. As such, it was primarily
intended to showcase hyperoptimization on large topology optimization problems. This is why it was built off
[TopOpt](https://www.topopt.mek.dtu.dk/apps-and-software/large-scale-topology-optimization-code-using-petsc), which is arguably the
highest performing large-scale topology optimization solver to date. As such, there are use cases HypOptLib is good for, and those it
is probably not best suited to.

### What is HypOptLib for?

As it was built off the TopOpt Petsc code, there are two main use cases that HypOptLib excels at. The first is any derivative topology optimization code
from topopt; ie anything TopOpt can do, HypOptLib can apply hyperoptimization to it. Most problems can be adapted in Python, and some more advanced cases might
need a little more work in the C++ code. But in general, with few modification HypOptLib will be compatible with any artbitrary topology optimization problem that
is supported by TopOpt.

The second use case that HypOptLib is suitable for is optimization problems that use Petsc in general. The Hyperoptimization class was abstracted from the
original TopOpt code, with wrappers provided for filtering, Lagrange Multipliers, and sensitivities. As such, problems that require the high-performance, scalable
solvers that Petsc provides for their sensitivity calculations can simply implement derived classes for each wrapper and use rest of HypOptLib unchanged.

### What is it not for?

Any optimization problem that does not fall into the above two categories is probably not suited to HypOptLib. The nature of the library is that it is focused on
large, parallelized workloads written explicitly for Petsc. While the classes are abstracted, the underlying data structures are all Petsc data structures. The
hyperoptimization algorithm is, however, very lightweight. This means that generalized hyperoptimization codes in Python, C++, or really any language would serve
much better for a large portion of optimization problems that can benefit from hyperoptimization.

## Installation and setup

HypOptLib is designed to run in Linux environments. It has been tested mainly
on WSL and Ubuntu, but should in theory support any Linux distribution with the
appropriate requirements.

HypOptLib can be setup in two main ways. A conveniance script `macros.sh` is
provided which can automatically install any prerequisites along with other
functionality. Alternatively, the prerequisites and their recommended method of
installation are also provided here.

### Automatic Install

The automatic install is only tested to work on Ubuntu. However it is not perfect and is mostly intended
as a shortcut to get started, you may still need to put in some leg work! It is especially **not**
recommended to run this if you do not own the system!! Always contact system administrators before
installing any libraries, and follow the procedures for the particular system or cluster you are using.

1. Clone the HypOptLib repository from https://github.com/Aidan-Sheedy/HypOptLib.

2. Run the automatic install macro.

```bash
./macros.sh setup
```

3. Restart the console

### Manual Install

For manual install, clone the HypOptLib repository as above, then install each
of the following requirements. You may have to run

```bash
sudo apt update
```

before installing dependencies.

| Dependency  | Command (Ubuntu)                                 |
|:-----------:|:-------------------------------------------------|
| cmake       | sudo apt install cmake                           |
| make        | sudo apt install make                            |
| mpi         | sudo apt install mpich                           |
| python pip3 | sudo apt install python3-pip                     |
| Pybind11    | pip3 install "pybind11[global]"                  |
| HDF5        | sudo apt install libhdf5-serial-dev              |
| BLAS        | sudo apt-get install libblas-dev liblapack-dev   |

Finally install PETSc, following the instructions here: https://petsc.org/release/install/download/.

When configuring PETSc, you must add the `--download-hdf5` tag:

```bash
./configure --download-hdf5
```

Some operating systems or environments may also need additional configure flags. If the compilation
fails, check the [PETSc install reference](https://petsc.org/release/install/install/) for help.
