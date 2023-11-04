
=====================================
HypOptLib documentation
=====================================

HypOptLib is a Petsc based semi-parallel implementation of hyperoptimization with a Python wrapper. The package is a
fork of the topology optimization framework
`TopOpt <https://www.topopt.mek.dtu.dk/apps-and-software/large-scale-topology-optimization-code-using-petsc>`_ by Aage
et. al., along with some borrowed elements from the derivative `TopOptLib <https://doi.org/10.1007/s00158-021-03018-7>`_ by Smit et. al.

Credit
========================
| **Hyperoptimization Algorithm** Hazhir Aliahmadi
| **Code** Aidan Sheedy
| **Advisor** Greg van Anders

Borrowed Material
========================
| **TopOpt** `Niels Aage, Erik Andreassen, Boyan Stefanov Lazarov <https://www.topopt.mek.dtu.dk/apps-and-software/large-scale-topology-optimization-code-using-petsc>`_
| **TopOptLib** `Thijs Smit, Niels Aage, Stephen J. Ferguson, Benedikt Halgason <https://doi.org/10.1007/s00158-021-03018-7>`_

Note
========================
This documentation will *not* cover the physics, derivation, or overview of the hyperoptimization algorithm.


.. toctree::
   :maxdepth: 2

   installation

.. toctree::
   :caption: Library Reference
   :maxdepth: 4
      
   api/library_root.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
