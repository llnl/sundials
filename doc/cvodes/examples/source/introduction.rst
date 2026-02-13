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

.. _CVODES.Examples.Intro:

Introduction
============

This document is the companion for the CVODES example programs.  It provides
explanations and sample outputs for the example programs that are distributed
with CVODES, and is intended to be used together with the CVODES User Guide.

Example types
-------------

The CVODES distribution contains examples of the following types:

- Serial and parallel examples of initial-value-problem (IVP) integration.
- Serial and parallel examples of forward sensitivity analysis (FSA).
- Serial and parallel examples of adjoint sensitivity analysis (ASA).
- Examples using OpenMP and Fortran interface modules.

Example names
-------------

Notes: example program names follow the pattern ``[slv][PbName]_[SA]_[ls]_[prec]_[p]``
where ``[slv]`` is the solver prefix (``cvs`` for CVODES examples), ``[SA]``
indicates sensitivity analysis (``FSA``, ``ASAi``, ``ASAp``), ``[ls]`` is the
linear-solver designation, ``[prec]`` is a preconditioner tag (if used), and
``[p]`` indicates the parallel (MPI) versions.

Below are the example program names as they appear in the distribution (grouped).

IVP (serial and parallel)
-------------------------

- ``cvsRoberts_dns``, ``cvsRoberts_dnsL``, ``cvsRoberts_dns_uw``
- ``cvsRoberts_klu``, ``cvsRoberts_sps``
- ``cvsAdvDiff_bnd``, ``cvsAdvDiff_bndL``
- ``cvsDiurnal_kry``, ``cvsDiurnal_kry_bp``
- ``cvsDirectDemo_ls``, ``cvsKrylovDemo_ls``, ``cvsKrylovDemo_prec``
- ``cvsParticle_dns``, ``cvsPendulum_dns``, ``cvsAnalytic_mels``
- ``cvs_analytic_fp_f2003``

FSA (forward sensitivity analysis)
----------------------------------

- ``cvsRoberts_FSA_dns``, ``cvsRoberts_FSA_dns_constraints``
- ``cvsRoberts_FSA_klu``, ``cvsRoberts_FSA_sps``
- ``cvsAdvDiff_FSA_non``, ``cvsDiurnal_FSA_kry``
- ``cvsRoberts_FSA_dns_Switch``, ``cvsAdvDiff_FSA_non_f2003``

ASA (adjoint sensitivity analysis)
----------------------------------

- ``cvsRoberts_ASAi_dns``, ``cvsRoberts_ASAi_dns_constraints``
- ``cvsRoberts_ASAi_klu``, ``cvsRoberts_ASAi_sps``
- ``cvsAdvDiff_ASAi_bnd``, ``cvsFoodWeb_ASAi_kry``
- ``cvsFoodWeb_ASAp_kry``, ``cvsHessian_ASA_FSA``

Fortran examples
----------------

A small set of examples are also provided using the SUNDIALS Fortran interface
modules; these are located in the ``examples/cvodes/F2003_serial`` and
``examples/cvodes/F2003_parallel`` directories and include, for example,
``cvs_analytic_fp_f2003`` and ``cvsAdvDiff_FSA_non_f2003``.

How to use these pages
----------------------

The following pages contain detailed descriptions of a subset of the examples,
sample output, and (where appropriate) figures illustrating the results.  Use
these examples as templates when integrating CVODES into your own application.
