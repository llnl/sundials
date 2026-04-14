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


.. _serial_f77:

====================================
Serial Fortran 77 example problems
====================================



.. _fark_diurnal_kry_bp:

fark_diurnal_kry_bp
===================================================

This problem is an ARKODE clone of the CVODE problem,
``fcv_diurnal_kry_bp``.  As described in :cite:p:`cvode_ug`, this problem
models a two-species diurnal kinetics advection-diffusion PDE system
in two spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} &=
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
problem as one involving :math:`2*Mx*My` coupled ODEs. In this
problem, we use a relatively coarse uniform mesh with
:math:`Mx=My=10`.

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver module, and the ARKSPILS interface.

The left preconditioner used is a banded matrix, constructed using
the ARKBP module.  The banded preconditioner matrix is generated using
difference quotients, with half-bandwidths ``mu = ml = 2``.

Performance data and sampled solution values are printed at
selected output times, and all performance counters are printed
on completion.






.. _fark_roberts_dnsL:

fark_roberts_dnsL
===================================================

This problem is an ARKODE clone of the CVODE problem,
``fcv_roberts_dnsL``.  As described in :cite:p:`cvode_ug`, this problem models
the kinetics of a three-species autocatalytic reaction.  This is an
ODE system with 3 components, :math:`Y = [y_1,\, y_2,\, y_3]^T`,
satisfying the equations,

.. math::

   \frac{d y_1}{dt} &= -0.04 y_1 + 10^4 y_2 y_3, \\
   \frac{d y_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3\cdot10^7 y_2^2, \\
   \frac{d y_3}{dt} &= 3\cdot10^7 y_2^2.

We integrate over the interval :math:`0\le t\le 4\cdot10^{10}`, with initial
conditions  :math:`Y(0) = [1,\, 0,\, 0]^T`.

Additionally, we supply the following two root-finding equations:

.. math::

   g_1(u) = u - 10^{-4}, \\
   g_2(w) = w - 10^{-2}.

While these are not inherently difficult nonlinear equations, they
easily serve the purpose of determining the times at which our
solutions attain desired target values.


Numerical method
----------------

This program solves the problem with a DIRK method, using a Newton
iteration with the SUNLINSOL_LAPACKDENSE linear solver module and
ARKDLS interface.

As with the :ref:`ark_robertson_root` problem, we enable ARKODE's
rootfinding module to find the times at which either :math:`u=10^{-4}`
or :math:`w=10^{-2}`.

Performance data and solution values are printed at
selected output times, along with additional output at rootfinding
events.  All performance counters are printed on completion.
