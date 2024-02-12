=======================================
Macros
=======================================

The macros.sh file is provided to simplify often-used sequences of commands, such
as building HypOptLib, documentation generation, etc. The following macros are currently
supported, where each macro is specified by calling::

    macros.sh <macro>

Supported Macros
=======================================

.. admonition:: -h

    Prints help string.

.. admonition:: build_docs

    Builds the documentation html. Requires the following packages:

    +------------+--------------------------------------------------+
    | Dependency | Command (Ubuntu)                                 |
    +============+==================================================+
    | Doxygen    | sudo apt-get install doxygen                     |
    +------------+--------------------------------------------------+
    | Sphinx     | sudo apt-get install python3-sphinx              |
    +------------+--------------------------------------------------+
    | Breathe    | Pip3 install breathe                             |
    +------------+--------------------------------------------------+
    | Sphinx RTD | pip3 install sphinx-rtd-theme                    |
    +------------+--------------------------------------------------+

.. admonition:: build

    Builds the HypOptLib package. The following options are supported:

        * **all** [default] will cmake and make commands
        * **clean** will cleanup all build files, cmake and make, and build them fresh.

.. admonition:: setup

    Sets up the build environment. This can be done manually as well. This does the following:

        1. Installs HypOptLib dependencies
        2. Installs PetSc dependencies
        3. Installs PetSc.
