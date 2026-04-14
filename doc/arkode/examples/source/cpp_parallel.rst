..
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _parallel_cpp:

====================================
Parallel C++ example problems
====================================



.. _ark_heat2D:

ark_heat2D
======================================================================

ARKODE provides one parallel C++ example problem, that extends our
previous :ref:`ark_heat1D` test to now simulate a two-dimensional heat
equation,

.. math::

   \frac{\partial u}{\partial t} = k_x \frac{\partial^2 u}{\partial x^2}
                                 + k_y \frac{\partial^2 u}{\partial y^2} + h,

for :math:`t \in [0, 0.3]`, and :math:`(x,y) \in [0, 1]^2`, with initial
condition :math:`u(0,x,y) = 0`, stationary boundary conditions,

.. math::

   \frac{\partial u}{\partial t}(t,0,y) = \frac{\partial u}{\partial t}(t,1,y) =
   \frac{\partial u}{\partial t}(t,x,0) = \frac{\partial u}{\partial t}(t,x,1) = 0,

and a periodic heat source,

.. math::

   h(x,y) = \sin(\pi x) \sin(2\pi y).

Under these conditions, the problem has an analytical solution of the
form

.. math::

   u(t,x,y) = \frac{1 - e^{-(k_x+4k_y)\pi^2 t}}{(k_x+4k_y)\pi^2} \sin(\pi x) sin(2\pi y).


Numerical method
----------------

The spatial derivatives are computed using second-order centered
differences, with the data distributed over :math:`nx\times ny` points
on a uniform spatial grid.

The problem is set up to use spatial grid parameters :math:`nx=60` and
:math:`ny=120`, with heat conductivity parameters :math:`k_x=0.5` and
:math:`k_y=0.75`.  The problem is run using scalar relative and
absolute solver tolerances of :math:`rtol=10^{-5}` and
:math:`atol=10^{-10}`.

As with the 1D version, this program solves the problem with a DIRK
method, that itself uses a Newton iteration and SUNLINSOL_PCG
iterative linear solver through the ARKSPILS interface.  However,
unlike the previous example, here the PCG solver is preconditioned
using a single Jacobi iteration, and uses ARKSPILS' built-in
finite-difference Jacobian-vector product routine. Additionally, this
problem uses MPI and the NVECTOR_PARALLEL module for parallelization.




Solutions
---------

Top row: 2D heat PDE solution snapshots, the left is at time :math:`t=0`,
center is at time :math:`t=0.03`, right is at time :math:`t=0.3`.
Bottom row, absolute error in these solutions.  Note that the relative
error in these results is on the order :math:`10^{-4}`, corresponding
to the spatial accuracy of the relatively coarse finite-difference
mesh.  All plots are created using the supplied Python script,
``plot_heat2D.py``.


.. image:: figs/plot-ark_heat2d_1.png
   :width: 30 %
.. image:: figs/plot-ark_heat2d_2.png
   :width: 30 %
.. image:: figs/plot-ark_heat2d_3.png
   :width: 30 %

.. image:: figs/plot-ark_heat2d_err_1.png
   :width: 30 %
.. image:: figs/plot-ark_heat2d_err_2.png
   :width: 30 %
.. image:: figs/plot-ark_heat2d_err_3.png
   :width: 30 %
