
========================
Getting Started
========================

HypOptLib is designed to run in Linux environments. It has been tested mainly
on WSL and Ubuntu, but should in theory support any Linux distribution with the
appropriate requirements.

Installation and Setup
========================

HypOptLib can be setup in two main ways. A conveniance script `macros.sh` is
provided which can automatically install any prerequisites along with other
functionality. Alternatively, the prerequisites and their recommended method of
installation are also provided here. For full documentation of the macros script
see this link **TODO - include link.**

Automatic Install
------------------------

1. Clone the HypOptLib repository from https://github.com/Aidan-Sheedy/Hyperoptimization_using_Petsc.git.

2. Run the automatic install macro:

    .. code-block:: bash

        ./macros.sh setup

Manual Install
------------------------

For manual install, clone the HypOptLib repository as above, then install each
of the following requirements:

+--------+---------------------------------+----------------+
| Symbol |           Description           | Units or Value |
+--------+---------------------------------+----------------+
|    1   | Energy of the incident Particle |        1       |
+--------+---------------------------------+----------------+
|    1   |        Mass per unit area       |        1       |
+--------+---------------------------------+----------------+
|    1   |                1                |        1       |
+--------+---------------------------------+----------------+
|    1   |                1                |        1       |
+--------+---------------------------------+----------------+
|    1   |                1                |        1       |
+--------+---------------------------------+----------------+


Building and Running
========================

They macros script again provides automatic building, but this can still be done
manually.
