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


.. _parallel_f77:

====================================
Parallel Fortran 77 example problems
====================================



.. _fark_diag_kry_bbd_p:

fark_diag_kry_bbd_p
===================================================

This problem is an ARKODE clone of the CVODE problem,
``fcv_diag_kry_bbd_p``.  As described in :cite:p:`cvode_ug`, this problem
models a stiff, linear, diagonal ODE system,

.. math::

   \frac{\partial y_i}{\partial t} &= -\alpha i y_i, \quad i=1,\ldots N.


Here :math:`\alpha=10` and :math:`N=10 N_P`, where :math:`N_P` is the
number of MPI tasks used for the problem.  The problem has initial
conditions :math:`y_i=1` and evolves for the time interval :math:`t\in
[0,1]`.




Numerical method
----------------

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned SUNLINSOL_SPGMR iterative linear
solver module and ARKSPILS interface.

A diagonal preconditioner matrix is used, formed automatically through
difference quotients within the ARKBBDPRE module.  Since ARKBBDPRE is
developed for use of a block-banded preconditioner, in this solver
each block is set to have half-bandwidths ``mudq = mldq = 0`` to
retain only the diagonal portion.

Two runs are made for this problem, first with left and then with
right preconditioning (``IPRE`` is first set to 1 and then to 2).

Performance data is printed at selected output times, and maximum
errors and final performance counters are printed on completion.







.. _fark_diag_non_p:

fark_diag_non_p
===================================================

This problem is an ARKODE clone of the CVODE problem,
``fcv_diag_non_p``.  As described in :cite:p:`cvode_ug`, this problem models a
nonstiff, linear, diagonal ODE system,

.. math::

   \frac{\partial y_i}{\partial t} &= -\alpha i y_i, \quad i=1,\ldots N.


Here :math:`\alpha=\frac{10}{N}` and :math:`N=10 N_P`, where :math:`N_P` is the
number of MPI tasks used for the problem.  The problem has initial
conditions :math:`y_i=1` and evolves for the time interval :math:`t\in [0,1]`.




Numerical method
----------------

This program solves the problem with an ERK method, and hence does not
require either a nonlinear or linear solver for integration.

Performance data is printed at selected output times, and maximum
errors and final performance counters are printed on completion.
