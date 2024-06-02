========================
Background
========================

Systematically understanding the solution space of large-scale design problems
is crucial for answering a range of design questions, including sensitivity,
performance, etc. These questions are particularly pressing for design problems
where the answer depends on a large number of design variables, e.g.,
high-resolution problems in topology optimization. HypOptLib is a C++/Python code
that implements novel filtration approaches to explore the solution space of
large-scale design problems. Rather than merely seeking an optimized solution,
HypOptLib yields detailed insight into the structure of the solution space, a
capability often overlooked in conventional optimization libraries. HypOptLib
exploits a mapping of design variables in conventional problems in topology
optimization to the dynamics of pseudothermal systems of “particles.” We
implement design variable dynamics by adapting techniques from molecular
dynamics. Using compliance minimization as a demonstration problem, we
demonstrate that HypOptLib supplements problems typically studied via topology
optimization with critical insights into the reliability and sensitivity of
design solutions. The problem-independent algorithm makes HypOptLib adaptable to
a wide range of large-scale design problems in topology optimization and beyond
where understanding. 

Methodology Overview
------------------------
Topology optimization concerns optimizing material layout within a given space,
aimed at meeting specific performance criteria under constraints such as material
volume. Central to our approach is the conceptualization of the design field
:math:`x`, which allows for the expression of design configurations via the number
of configurations :math:`\Omega(C,V)`, facilitated through the application of the
Pareto-Laplace filter [1]:

.. math::
  Z(\beta) = \int_{C_\text{min}(V)}^{\infty}dC e^{-\beta C} \Omega(C,V),

where :math:`\beta` represents the Laplace variable and :math:`C_\text{min}` denotes the
minimal compliance for a given material volume :math:`V`. :math:`Z(\beta)` thus acts as a
generating function, encoding the solution space's geometry through a weighted
summation over possible design realizations.

Direct evaluation of :math:`Z(\beta)` is impractical due to unknown variables such as
:math:`C_\text{min}(V)` and :math:`\Omega(C,V)`. To navigate this, we adopt practices from
statistical physics, leveraging underlying degrees of freedom within design
fields to infer these quantities. This methodology enables a finite element
approximation, simplifying :math:`Z(\beta)` as follows:

.. math::
    Z(\beta) = \sum_{\{x\}} e^{-\beta C} \delta\left(\sum_{e}x_e - V\right),

where the sum encompasses all material distributions :math:`\{x_e\}` constrained by
material volume through the Dirac delta function.

This redefined :math:`Z(\beta)` serves to geometrize the compliance minimization
solution space, with :math:`\Omega(C,V)` representing the material distribution
patterns meeting specific compliance and volume criteria. As :math:`\beta` varies,
we transition from emphasizing configurations near minimal compliance to treating
all configurations equally, providing a comprehensive understanding of the
solution space dynamics [2].

Numerical Implementation and Molecular Dynamics
_______________________________________________
The numerical implementation relies on identifying :math:`Z(\beta)` as a partition
function, akin to those in molecular dynamics, treating compliance as potential
energy and material cells as "particles" within the design domain. This analogy
extends to computing :math:`Z(\beta)` by introducing auxiliary variables representing
particle momentum, integrating kinetic energy considerations to maintain dynamic
configurations within the solution space [2].

Implementing a Nosé-Hoover chain thermostat allows for simulating these dynamics
under constant temperature, ensuring that the generated configurations
appropriately reflect the integral transform of the compliance minimization
problem. This approach, detailed further in our provided open-source
implementations, leverages conventional molecular dynamics techniques proven
effective in large-scale systems.

The methodology concludes with the formulation of equations of motion, employing
a symplectic integrator for stable numerical simulation. This framework ensures
that generated configurations respect both compliance and material volume
constraints, offering a robust and scalable solution to exploring the solution
space in topology optimization.

References
---------------------
[1] Aliahmadi, H., Perez, R. and van Anders, G. (2024) Transforming Design Spaces
Using Pareto-Laplace Filters [Preprint]. doi:arXiv:2403.00631.

[2] Aliahmadi, H., Sheedy, A., Perez, R. and van Anders, G. (2024) Hyperoptimization insight
for computational morphogenesis [Preprint]. [Unpublished]
