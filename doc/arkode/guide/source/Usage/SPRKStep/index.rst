.. ----------------------------------------------------------------
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

.. _ARKODE.Usage.SPRKStep:

==========================================
Using the SPRKStep time-stepping module
==========================================

This section is concerned with the use of the SPRKStep time-stepping
module for the solution of initial value problems (IVPs) in a C or C++
language setting.  Usage of SPRKStep follows that of the rest of ARKODE,
and so in this section we primarily focus on those usage aspects that
are specific to SPRKStep.

We note that of the ARKODE example programs located in the source code
``examples/arkode`` folder, the following demonstrate ``SPRKStep`` usage:

* ``examples/arkode/C_serial/ark_harmonic_symplectic.c``
* ``examples/arkode/C_serial/ark_damped_harmonic_symplectic.c``, and
* ``examples/arkode/C_serial/ark_kepler.c``

.. toctree::
   :maxdepth: 1

   User_callable
