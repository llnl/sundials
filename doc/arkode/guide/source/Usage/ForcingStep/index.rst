.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
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

.. _ARKODE.Usage.ForcingStep:

==========================================
Using the ForcingStep time-stepping module
==========================================

This section is concerned with the use of the ForcingStep time-stepping module
for the solution of initial value problems (IVPs) in a C or C++ language
setting.  Usage of ForcingStep follows that of the rest of ARKODE, and so in
this section we primarily focus on those usage aspects that are specific to
ForcingStep. A skeleton of a program using ForcingStep follows essentially the
same structure as SplittingStep
(see :numref:`ARKODE.Usage.SplittingStep.Skeleton`).

.. toctree::
   :maxdepth: 1

   User_callable
