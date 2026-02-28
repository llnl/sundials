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

.. _openmp_c:

====================================
OpenMP C example problems
====================================




.. _ark_brusselator1D_omp:

ark_brusselator1D_omp
============================================

This problem is mathematically identical to the one-dimensional
reaction-diffusion brusselator model, :ref:`ark_brusselator1D`.  As
before, we investigate a time-dependent system of partial differential
equations with 3 components, :math:`Y = [u,\, v,\, w]^T` that satisfy
the equations,

.. math::

   \frac{\partial u}{\partial t} &= d_u \frac{\partial^2 u}{\partial
      x^2} + a - (w+1) u + v u^2, \\
   \frac{\partial v}{\partial t} &= d_v \frac{\partial^2 v}{\partial
      x^2} + w u - v u^2, \\
   \frac{\partial w}{\partial t} &= d_w \frac{\partial^2 w}{\partial
      x^2} + \frac{b-w}{\varepsilon} - w u.

We integrate for :math:`t \in [0, 10]`, and :math:`x \in [0, 1]`, with
initial conditions

.. math::

   u(0,x) &=  a + \frac{1}{10} \sin(\pi x),\\
   v(0,x) &= \frac{b}{a} + \frac{1}{10}\sin(\pi x),\\
   w(0,x) &=  b + \frac{1}{10}\sin(\pi x),

and with stationary boundary conditions, i.e.

.. math::

   \frac{\partial u}{\partial t}(t,0) &= \frac{\partial u}{\partial t}(t,1) = 0,\\
   \frac{\partial v}{\partial t}(t,0) &= \frac{\partial v}{\partial t}(t,1) = 0,\\
   \frac{\partial w}{\partial t}(t,0) &= \frac{\partial w}{\partial t}(t,1) = 0.



Numerical method
----------------

The numerical method is identical to the previous implementation,
except that we now use SUNDIALS' OpenMP-enabled vector kernel module,
NVECTOR_OPENMP, and have similarly threaded the supplied right-hand
side and banded Jacobian construction functions.
