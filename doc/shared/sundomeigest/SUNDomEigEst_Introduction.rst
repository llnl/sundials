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

.. _SUNDomEigEst.Introduction:

Introduction to Dominant Eigenvalue Estimators
==============================================

For problems that require the dominant eigenvalue of a matrix (i.e., the Jacobian),
the SUNDIALS packages operate using generic dominant eigenvalue estimator modules
defined through the :c:type:`SUNDomEigEstimator` class.
This allows SUNDIALS packages to utilize any valid :c:type:`SUNDomEigEstimator`
implementation that provides a set of required functions.  These
functions can be divided into three categories.  The first are the core
estimator functions.  The second group consists of "set" routines
to supply the dominant eigenvalue estimator object with functions provided by the
SUNDIALS package, or for modification of estimator parameters.  The last
group consists of "get" routines for retrieving artifacts (statistics,
residual, etc.) from the estimator.  All of these functions
are defined in the header file ``sundials/sundials_domeigestimator.h``.

The implementations provided with SUNDIALS work in coordination
with the SUNDIALS :c:type:`N_Vector` modules to provide a set of compatible data
structures for the estimator.
Moreover, advanced users can provide a customized :c:type:`SUNDomEigEstimator`
implementation to any SUNDIALS package, particularly in cases where they
provide their own :c:type:`N_Vector`.

While Krylov-based estimators preset the number of Krylov subspace
dimensions, resulting in a tolerance-free estimation, SUNDIALS requires
that iterative estimators stop when the residual meets a prescribed
tolerance, :math:`\tau`,

.. math::
  :name: pi_rel_tol

  \frac{\left|\lambda_k - \lambda_{k-1}\right|}{\left|\lambda_k \right|} < \tau.

For users interested in providing their own :c:func:`SUNDomEigEstimator`, the
following section presents the :c:type:`SUNDomEigEstimator` class and its implementation
beginning with the definition of :c:type:`SUNDomEigEstimator` functions in
:numref:`SUNDomEigEst.CoreFn` -- :numref:`SUNDomEigEst.GetFn`. This is followed by
the definition of functions supplied to an estimator implementation in
:numref:`SUNDomEigEst.SUNSuppliedFn`. The :c:type:`SUNDomEigEstimator` type is defined
:numref:`SUNDomEigEst.Generic`. The section that then follows describes
the :c:type:`SUNDomEigEstimator` functions required by this SUNDIALS package, and provides
additional package specific details. Then the remaining sections of this
chapter present the :c:type:`SUNDomEigEstimator` modules provided with SUNDIALS.
