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


.. _serial_cpp:

====================================
Serial C++ example problems
====================================


.. _ark_analytic_sys:

ark_analytic_sys
===============================================

This example demonstrates the use of ARKODE's fully implicit solver on
a stiff ODE system that has a simple analytical solution.  The problem
is that of a linear ODE system,

.. math::

   \frac{dy}{dt} = Ay

where :math:`A = V D V^{-1}`.  In this example, we use

.. math::

   V = \left[\begin{array}{rrr} 1 & -1 & 1\\ -1 & 2 & 1\\ 0 & -1 & 2
       \end{array}\right], \qquad
   V^{-1} = \frac14 \left[\begin{array}{rrr} 5 & 1 & -3\\ 2 & 2 & -2\\
       1 & 1 & 1 \end{array}\right], \qquad
   D = \left[\begin{array}{rrr} -1/2 & 0 & 0\\ 0 & -1/10 & 0\\ 0 & 0 &
       \lambda \end{array}\right].

where :math:`\lambda` is a large negative number. The analytical
solution to this problem may be computed using the matrix exponential,

.. math::

   Y(t) = V e^{Dt} V^{-1} Y(0).

We evolve the problem for :math:`t` in the interval :math:`\left[0,\,
\frac{1}{20}\right]`, with initial condition :math:`Y(0) = \left[1,\,
1,\, 1\right]^T`.


Numerical method
----------------

The stiffness of the problem is directly proportional to the
value of :math:`\lambda`.  The value of :math:`\lambda` should be
negative to result in a well-posed ODE; for values with magnitude
larger than 100 the problem becomes quite stiff.

Here, we choose :math:`\lambda = -100`, along with scalar relative and
absolute tolerances of :math:`rtol=10^{-6}` and :math:`atol=10^{-10}`,
respectively.

This program solves the problem with the DIRK method,
Newton iteration with the SUNMATRIX_DENSE matrix module and
accompanying SUNLINSOL_DENSE linear solver module, ARKDLS direct
linear solver interface, and a user-supplied dense Jacobian
routine.  Output is printed every 0.005 units of time (10 total).
Run statistics (optional outputs) are printed at the end.



Solutions
---------

This problem is included both as a simple example to test systems of
ODE within ARKODE on a problem having an analytical solution,
:math:`Y(t) = V e^{Dt} V^{-1} Y(0)`.  As seen in the plots below, the
computed solution tracks the analytical solution quite well (left),
and results in errors with exactly the magnitude as specified by the
requested error tolerances (right).

.. image:: figs/plot-ark_analytic_sys.png
   :width: 45 %
.. image:: figs/plot-ark_analytic_sys_error.png
   :width: 45 %
