..
   -----------------------------------------------------------------------------
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
   -----------------------------------------------------------------------------

.. _CVODES.Examples.Adjoint:

Adjoint sensitivity analysis examples
=====================================

A serial dense example: cvsRoberts_ASAi_dns
-------------------------------------------

This example demonstrates the use of CVODES adjoint sensitivity analysis to
compute gradients of a time-integrated functional with respect to model
parameters without explicitly computing solution sensitivities.

.. literalinclude:: ../../../../examples/cvodes/serial/cvsRoberts_ASAi_dns.out
   :language: none

A parallel adjoint example: cvsAdvDiff_ASAp_non_p
-------------------------------------------------

This example shows a parallel adjoint calculation for a 1-D advection-diffusion
problem and the computation of time-integrated gradients via appended
quadrature equations.

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsAdvDiff_ASAp_non_p.out
   :language: none

A parallel adjoint example with CVBBDPRE: cvsAtmDisp_ASAi_kry_bbd_p
-------------------------------------------------------------------

This example models an atmospheric dispersion problem and uses the CVBBDPRE
preconditioner with Krylov solvers in both forward and backward phases.

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsAtmDisp_ASAi_kry_bbd_p.out
   :language: none
