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

.. _CVODES.Examples.C:

C Example Problems
==================

A serial dense example: cvsRoberts_dns
--------------------------------------

This example solves a simple three-species chemical kinetics problem and
illustrates the use of the BDF method with Newton iteration and a dense
linear solver.  See the original LaTeX source for more details.

.. literalinclude:: ../../../../examples/cvodes/serial/cvsRoberts_dns.out
   :language: none

A serial dense forward-sensitivity example: cvsRoberts_FSA_dns
--------------------------------------------------------------

This example modifies the chemical kinetics problem to compute forward
sensitivities with respect to reaction rate constants using CVODES forward
sensitivity capabilities.

.. literalinclude:: ../../../../examples/cvodes/serial/cvsRoberts_FSA_dns_-sensi_sim_t.out
   :language: none

A parallel example: cvsDiurnal_FSA_kry_p
----------------------------------------

This parallel example solves a 2-D advection-diffusion PDE and demonstrates
use of the MPI NVECTOR and Krylov linear solvers together with a block
preconditioner.

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsDiurnal_FSA_kry_p_-sensi_sim_t.out
   :language: none


.. note::

   These reST pages are mechanically generated from the original LaTeX
   example documentation.  Some LaTeX constructs and custom macros were
   converted heuristically and may require manual edits for perfect Sphinx
   rendering.  See the conversion notes in the Pull Request for details.
