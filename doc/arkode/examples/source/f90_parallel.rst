..
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
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


.. _parallel_f90:

====================================
Parallel Fortran 90 example problems
====================================


.. _fark_heat2D:

fark_heat2D
===================================================

This test problem is a Fortran-90 version of the same two-dimensional
heat equation problem as in C++, :ref:`ark_heat2D`.  This models a
simple two-dimenaional heat equation,

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

The spatial grid is set to :math:`nx=60` and :math:`ny=120`.  The heat
conductivity parameters are :math:`k_x=0.5` and :math:`k_y=0.75`.

As with the C++ version, this program solves the problem with a DIRK
method, that itself uses a Newton iteration and SUNLINSOL_PCG
iterative linear solver module through the ARKSPILS interface.  Also
like the C++ version, the PCG solver is preconditioned using a single
Jacobi iteration, and uses ARKSPILS' finite-difference Jacobian-vector
product approximation routine for the PCG polver.  Additionally, this
problem uses MPI and the Fortran interface to the NVECTOR_PARALLEL
module for parallelization.
