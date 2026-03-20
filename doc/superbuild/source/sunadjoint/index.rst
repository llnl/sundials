..
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

.. _SUNAdjoint:

############################
Adjoint Sensitivity Analysis
############################

This section presents the :c:type:`SUNAdjointStepper` and
:c:type:`SUNAdjointCheckpointScheme` classes that provide a common interface for
adjoint sensitivity analysis (ASA) capabilities. Currently it supports :ref:`the
ASA capabilities in ARKODE <ARKODE.Mathematics.ASA>`, while the ASA capabilities
in :ref:`CVODES <CVODES.Mathematics.ASA>` and :ref:`IDAS <IDAS.Mathematics.ASA>`
must be used directly.

.. toctree::
   :maxdepth: 1

   SUNAdjoint_links.rst
