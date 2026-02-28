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

.. _KINSOL.Examples.Intro:

Introduction
============

KINSOL includes examples of many types, to illustrate the use of various
nonlinear and linear solver options and vector implementations.

With the exception of "demo"-type example files, the names of the examples are
generally of the form ``[slv][PbName]_[strat]_[ls]_[prec]_[p]``, where:

- ``[slv]`` identifies the solver (in this case ``kin``).
- ``[PbName]`` identifies the problem.
- ``[strat]`` identifies the strategy (absent if "none" or "linesearch").
- ``[ls]`` identifies the linear solver module used.
- ``[prec]`` indicates the KINSOL preconditioner module used (only if
  applicable, for examples using a Krylov linear solver and the ``KINBBDPRE``
  module, this will be ``bbd``).
- ``[p]`` indicates an example using the MPI parallel vector.

The following lists summarize all the examples distributed with KINSOL.  In the
subsequent sections, we give detailed descriptions of some (but not all) of
these examples. We also give our output files for each of these examples, but
users should be cautioned that their results may differ slightly from
these. Differences in solution values may differ within the tolerances, and
differences in cumulative counters, such as numbers of Newton iterations, may
differ from one machine environment to another by as much as 10% to 20%.

.. note::

   The examples in the KINSOL distribution are written in such a way as to
   compile and run for any combination of configuration options used during the
   installation of SUNDIALS. As a consequence, they contain portions of code
   that will not be typically present in a user program. For example, programs
   may make use of the variables ``SUNDIALS_EXTENDED_PRECISION`` and
   ``SUNDIALS_DOUBLE_PRECISION`` to test if the solver libraries were built in
   extended or double precision, and use the appropriate conversion specifiers
   in ``printf`` functions.

Serial examples
---------------

Supplied in the ``examples/kinsol/serial`` directory are the following serial
examples (using the :ref:`serial vector <NVectors.NVSerial>`):

- ``kinRoberts_fp`` solves the backward Euler time step for a three-species
  chemical kinetics system, using the fixed point strategy.

- ``kinFerTron_dns`` solves the Ferraris-Tronconi problem.

  This program solves the problem with the :ref:`dense linear solver
  <SUNLinSol_Dense>` and uses different combinations of globalization and
  Jacobian update strategies with different initial guesses.

- ``kinFerTron_klu`` solves the same problem as in ``kinFerTron_dns``, but uses
  the :ref:`KLU sparse direct linear solver <SUNLinSol.KLU>`.

- ``kinRoboKin_dns`` solves a nonlinear system from robot kinematics.

  This program solves the problem with the :ref:`dense linear solver
  <SUNLinSol_Dense>` and a user-supplied Jacobian routine.

- ``kinRoboKin_slu`` is the same as ``kinRoboKin_dns`` but uses the
  :ref:`SuperLU_MT sparse direct linear solver <SUNLinSol.SuperLUMT>`.

- ``kinLaplace_bnd`` solves a simple 2D elliptic PDE on a unit square.

  This program solves the problem with the :ref:`banded direct linear solver
  <SUNLinSol_Band>`.

- ``kinLaplace_picard_bnd`` is the same as ``kinLaplace_bnd`` but uses the
  Picard strategy.

- ``kinLaplace_picard_kry`` is the same as ``kinLaplace_picard_bnd`` but uses
  the :ref:`GMRES iterative linear solver <SUNLinSol.SPGMR>`.

- ``kinFoodWeb_kry`` solves a food web model.

  This is a nonlinear system that arises from a system of partial differential
  equations describing a six-species food web population model, with
  predator-prey interaction and diffusion on the unit square in two
  dimensions. This program solves the problem with the the :ref:`GMRES iterative
  linear solver <SUNLinSol.SPGMR>` and a user-supplied preconditioner. The
  preconditioner is a block-diagonal matrix based on the partial derivatives of
  the interaction terms only.

- ``kinKrylovDemo_ls`` solves the same problem as ``kinFoodWeb_kry``, but with
  three Krylov linear solvers: :ref:`GMRES <SUNLinSol.SPGMR>`, :ref:`BiCGSTAB
  <SUNLinSol.SPBCGS>`, and :ref:`TFQMR <SUNLinSol.SPTFQMR>`.

- ``kinAnalytic_fp.c`` solves a small nonlinear system with known solution using
  the fixed-point iteration.

MPI examples
------------

Supplied in the ``examples/kinsol/parallel`` directory are the following
parallel examples (using the :ref:`MPI parallel vector <NVectors.NVParallel>`):

- ``kinFoodWeb_kry_p`` is a parallel implementation of ``kinFoodWeb_kry``.

- ``kinFoodWeb_kry_bbd_p`` solves the same problem as ``kinFoodWeb_kry_p``, with
  the :ref:`KINBBDPRE <KINSOL.Usage.CC.kin_bbdpre>` band-block-diagonal
  preconditioner.

OpenMP example
--------------

Supplied in the ``examples/kinsol/C_openmp`` directory is an example
``kinFoodweb_kry_omp``, which solves the same problem as ``kinFoodweb_kry``, but
using the :ref:`OpenMP vector <NVectors.OpenMP>`.

MPI + CUDA example
------------------

Supplied in the ``examples/kinsol/CUDA_mpi`` directory is an example
``kin_em_mpicuda``, which solves an expectation-maximization problem for mixture
densities, using the :ref:`MPI+X <NVectors.MPIPlusX>` and :ref:`CUDA
<NVectors.CUDA>` vectors.

C++ examples
------------

- The directory ``examples/kinsol/CXX_parallel`` contains two
  examples. ``kin_heat2D_nonlin_p`` solves a steady state 2D heat equation with
  an additional nonlinear term. This example is solved with fixed point
  iteration and illustrates the use of various orthogonalization methods with
  Anderson acceleration. ``kin_em_p`` solves the same problem as
  ``kin_em_mpicuda`` listed above.

- The directory ``examples/kinsol/CXX_parhyp`` contains two examples,
  ``kin_heat2D_nonlin_hypre_pfmg`` and ``kin_bratu2D_hyper_pfmg``. These use the
  *hypre* PFMG preconditioner, fixed point iteration, and Anderson acceleration.

Fortran examples
----------------

The following are examples in Fortran, using the SUNDIALS Fortran interface
modules, and are based on examples listed earlier:

- The directory ``examples/kinsol/F2003_serial`` contains four examples:
  ``kinDiagon_kry_f2003`` solves a simple diagonal test
  problem. ``kinLaplace_bnd_f2003`` solves the same problem as
  ``kinLaplace_bnd``. ``kinLaplace_picard_kry_f2003`` solves the same problem as
  ``kinLaplace_picard_kry``. ``kinRoboKin_dns_f2003`` solves the same problem as
  ``kinRoboKin_dns``.

- The directory ``examples/kinsol/F2003_parallel`` contains one example,
  ``kin_diagon_kry_f2003``, which solves the same problem as
  ``kinDiagon_kry_f2003``, but using MPI.
