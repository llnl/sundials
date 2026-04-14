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

.. _parallel_c:

====================================
Parallel C example problems
====================================



.. _ark_diurnal_kry_bbd_p:

ark_diurnal_kry_bbd_p
===================================================


This problem is an ARKODE clone of the CVODE problem,
``cv_diurnal_kry_bbd_p``.  As described in :cite:p:`cvode_ug`, this problem
models a two-species diurnal kinetics advection-diffusion PDE system
in two spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} =
     K_h \frac{\partial^2 c_i}{\partial x^2} +
     V \frac{\partial     c_i}{\partial x} +
     \frac{\partial}{\partial y}\left( K_v(y)
     \frac{\partial c_i}{\partial y}\right) +
     R_i(c_1,c_2,t),\quad i=1,2

where

.. math::

   R_1(c_1,c_2,t) &= -q_1*c_1*c_3 - q_2*c_1*c_2 + 2*q_3(t)*c_3 + q_4(t)*c_2, \\
   R_2(c_1,c_2,t) &=  q_1*c_1*c_3 - q_2*c_1*c_2 - q_4(t)*c_2, \\
   K_v(y) &= K_{v0} e^{y/5}.

Here :math:`K_h`, :math:`V`, :math:`K_{v0}`, :math:`q_1`, :math:`q_2`,
and :math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)`
vary diurnally.  The problem is posed on the square spatial domain
:math:`(x,y) \in [0,20]\times[30,50]`, with homogeneous Neumann
boundary conditions, and for time interval :math:`t\in [0,86400]` sec
(1 day).

We enforce the initial conditions

.. math::

   c_1(x,y) &=  10^6 \chi(x)\eta(y) \\
   c_2(x,y) &=  10^{12} \chi(x)\eta(y) \\
   \chi(x) &= 1 - \sqrt{\frac{x - 10}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}} \\
   \eta(y) &= 1 - \sqrt{\frac{y - 40}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}}.




Numerical method
----------------

We employ a method of lines approach, wherein we first
semi-discretize in space to convert the system of 2 PDEs into a larger
system of ODEs.  To this end, the spatial derivatives are computed
using second-order centered differences, with the data distributed
over :math:`Mx*My` points on a uniform spatial grid.  As a result, ARKODE
approaches the problem as one involving :math:`2*Mx*My` coupled ODEs.

The problem is decomposed in parallel into uniformly-sized subdomains,
with two subdomains in each direction (four in total), and where each
subdomain has five points in each direction (i.e. :math:`Mx=My=10`).

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver through the ARKSPILS interface.

The preconditioner matrix used is block-diagonal, with banded blocks,
constructed using the ARKBBDPRE module.  Each block is generated using
difference quotients, with half-bandwidths ``mudq = mldq = 10``, but
the retained banded blocks have half-bandwidths ``mukeep = mlkeep = 2``.
A copy of the approximate Jacobian is saved and conditionally reused
within the preconditioner routine.

Two runs are made for this problem, first with left and then with
right preconditioning.

Performance data and sampled solution values are printed at
selected output times, and all performance counters are printed
on completion.




.. _ark_diurnal_kry_p:

ark_diurnal_kry_p
===================================================

This problem is an ARKODE clone of the CVODE problem,
``cv_diurnal_kry_p``.  As described in :cite:p:`cvode_ug`, this test problem
models a two-species diurnal kinetics advection-diffusion PDE system
in two spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} =
     K_h \frac{\partial^2 c_i}{\partial x^2} +
     V \frac{\partial     c_i}{\partial x} +
     \frac{\partial}{\partial y}\left( K_v(y)
     \frac{\partial c_i}{\partial y}\right) +
     R_i(c_1,c_2,t),\quad i=1,2

where

.. math::

   R_1(c_1,c_2,t) &= -q_1*c_1*c_3 - q_2*c_1*c_2 + 2*q_3(t)*c_3 + q_4(t)*c_2, \\
   R_2(c_1,c_2,t) &=  q_1*c_1*c_3 - q_2*c_1*c_2 - q_4(t)*c_2, \\
   K_v(y) &= K_{v0} e^{y/5}.

Here :math:`K_h`, :math:`V`, :math:`K_{v0}`, :math:`q_1`, :math:`q_2`,
and :math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)`
vary diurnally.  The problem is posed on the square spatial domain
:math:`(x,y) \in [0,20]\times[30,50]`, with homogeneous Neumann
boundary conditions, and for time interval :math:`t\in [0,86400]` sec
(1 day).

We enforce the initial conditions

.. math::

   c^1(x,y) &=  10^6 \chi(x)\eta(y) \\
   c^2(x,y) &=  10^{12} \chi(x)\eta(y) \\
   \chi(x) &= 1 - \sqrt{\frac{x - 10}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}} \\
   \eta(y) &= 1 - \sqrt{\frac{y - 40}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}}.




Numerical method
----------------

We employ a method of lines approach, wherein we first semi-discretize
in space to convert the system of 2 PDEs into a larger system of ODEs.
To this end, the spatial derivatives are computed using second-order
centered differences, with the data distributed over :math:`Mx*My`
points on a uniform spatial grid.  As a result, ARKODE approaches the
problem as one involving :math:`2*Mx*My` coupled ODEs.

The problem is decomposed in parallel into uniformly-sized subdomains,
with two subdomains in each direction (four in total), and where each
subdomain has five points in each direction (i.e. :math:`Mx=My=10`).

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver, through the ARKSPILS interface.

The preconditioner matrix used is block-diagonal, with block-diagonal
portion of the Newton matrix used as a left preconditioner.  A copy of
the block-diagonal portion of the Jacobian is saved and conditionally
reused within the preconditioner routine.

Performance data and sampled solution values are printed at
selected output times, and all performance counters are printed
on completion.
