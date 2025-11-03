.. ----------------------------------------------------------------
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

********************
ARKODE Documentation
********************

This is the documentation for ARKODE, an adaptive step time
integration package for stiff, nonstiff and mixed stiff/nonstiff
systems of ordinary differential equations (ODEs) using Runge--Kutta
(i.e., one-step, multi-stage) methods.  The ARKODE solver is a
component of the `SUNDIALS
<https://computing.llnl.gov/projects/sundials>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computing.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, including user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, and for calculations in serial, using shared-memory
parallelism (e.g., via OpenMP, CUDA, Raja, Kokkos) or distributed-memory
parallelism (via MPI).  The default integration and solver options
should apply to most users, though control over nearly all internal
parameters and time adaptivity algorithms is enabled through optional
interface routines.

ARKODE is written in C, with C++ and Fortran interfaces.

ARKODE is developed by `Southern Methodist University <https://www.smu.edu>`_ and
`Lawrence Livermore National Security <https://www.llnl.gov>`_,
with support by the `US Department of Energy <http://www.doe.gov>`_,
`Office of Science <https://www.energy.gov/science/office-science>`_.


.. include:: Landing.rst

.. only:: html

   **Table of Contents**


.. toctree::
   :numbered:
   :maxdepth: 1

   Introduction
   Mathematics
   Organization
   sundials/index.rst
   Usage/index.rst
   ARKodeButcherTable
   ARKodeSPRKTable
   nvectors/index.rst
   sunmatrix/index.rst
   sunlinsol/index.rst
   sunnonlinsol/index.rst
   sundomeigest/index.rst
   sunadaptcontroller/index.rst
   sunstepper/index.rst
   sunadjoint/index.rst
   sunmemory/index.rst
   sundials/Install_link.rst
   Constants
   Butcher
   Fortran/index.rst
   History_link.rst
   Changelog_link.rst
   References
.. only:: html

   * :ref:`genindex`
