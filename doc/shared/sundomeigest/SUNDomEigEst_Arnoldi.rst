..
   Programmer(s): Mustafa Aggul @ SMU
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

.. _SUNDomEigEst.Arnoldi:

The SUNDomEigEstimator_Arnoldi Module
=====================================

.. versionadded:: 7.5.0

The SUNDomEigEstimator_Arnoldi implementation of the
:c:type:`SUNDomEigEstimator` class performs the Arnoldi Iteration method
:cite:p:`arnoldi51`; this is an iterative dominant eigenvalue estimator that is
designed to be compatible with any ``N_Vector`` implementation that supports a
minimal subset of operations (:c:func:`N_VClone()`, :c:func:`N_VDotProd()`,
:c:func:`N_VScale()`, and :c:func:`N_VDestroy()`).

Arnoldi iteration is particularly effective for large, sparse matrices where only
the dominant eigenvalue is needed.  It constructs an orthonormal basis of the Krylov
subspace

.. math::

   \mathcal{K}_m(A, \mathbf{v}) = \text{span}\{\mathbf{v}, A \mathbf{v}, A^2 \mathbf{v}, \dots, A^{m-1} \mathbf{v}\}

using the Gram-Schmidt process.  The matrix :math:`A` is projected onto this subspace
to form a small upper Hessenberg matrix :math:`H_m`.  The eigenvalues of :math:`H_m`
approximate some of the eigenvalues of :math:`A`; the dominant eigenvalue of :math:`A` is
well-approximated by the dominant eigenvalue of :math:`H_m`.

Arnoldi iteration works for matrices with both real and complex eigenvalues.  It supports
estimations with a user-specified fixed Krylov subspace dimension (at least 3).  While
the choice of dimension results in a prefixed amount of memory, it strictly
determines the quality of the estimate.  To improve the estimation accuracy, we have found that
preprocessing with a number of Power iterations is particularly useful.
This operation is free from any additional memory requirement and is further explained below.

The matrix :math:`A` is not required explicitly; only a routine that provides an
approximation of the matrix-vector product, :math:`Av`, is required.


.. _SUNDomEigEst.Arnoldi.Usage:

SUNDomEigEstimator_Arnoldi Usage
--------------------------------

To use SUNDomEigEstimator_Arnoldi include the header file
``sundomeigest/sundomeigest_arnoldi.h``, and link to the library
``libsundials_sundomeigestarnoldi``.

The module SUNDomEigEstimator_Arnoldi provides the following user-callable
routines:


.. c:function:: SUNDomEigEstimator SUNDomEigEstimator_Arnoldi(N_Vector q, int kry_dim, SUNContext sunctx);

   This constructor function creates and allocates memory for the Arnoldi
   iteration implementation of a :c:type:`SUNDomEigEstimator`.

   Consistency checks are performed to ensure the input vector is non-zero and
   supplies the necessary operations.

   :param q: the initial guess for the dominant eigenvector; this should not be
             a non-dominant eigenvector of the Jacobian.
   :param kry_dim: the dimension of the Krylov subspace (default 3). A value
                   :math:`\leq 2` will result in using default value. This
                   default is chosen to minimize the memory footprint.
   :param sunctx: the :c:type:`SUNContext` object.

   :returns: If successful, a :c:type:`SUNDomEigEstimator` otherwise ``NULL``.

   .. note::

      When used in a time-dependent context, the initial guess supplied to the
      constructor, ``q``, is used only in the first
      :c:func:`SUNDomEigEstimator_Estimate` call and is overwritten with the
      result of the most recent preprocessing iterations (see
      :c:func:`SUNDomEigEstimator_SetNumPreprocessIters`). As an initial guess
      too close to the dominant eigenvector may cause a breakdown in the
      Gram–Schmidt process within the Arnoldi iteration, users should account
      for this when setting the number of initial and subsequent preprocessing
      iterations (e.g., with LSRKStep see
      :c:func:`LSRKStepSetNumDomEigEstInitPreprocessIters` and
      :c:func:`LSRKStepSetNumDomEigEstPreprocessIters`).

      The initial guess can be reset with
      :c:func:`SUNDomEigEstimator_SetInitialGuess`.


.. _SUNDomEigEst.Arnoldi.Description:

SUNDomEigEstimator_Arnoldi Description
--------------------------------------

The SUNDomEigEstimator_Arnoldi module defines the *content* field of a
``SUNDomEigEstimator`` to be the following structure:

.. code-block:: c

   struct SUNDomEigEstimatorContent_Arnoldi_ {
     SUNATimesFn ATimes;
     void* ATdata;
     N_Vector* V;
     N_Vector q;
     int kry_dim;
     int num_warmups;
     long int num_iters;
     long int num_ATimes;
     sunrealtype* LAPACK_A;
     sunrealtype* LAPACK_wr;
     sunrealtype* LAPACK_wi;
     sunrealtype* LAPACK_work;
     snuindextype LAPACK_lwork;
     sunrealtype** LAPACK_arr;
     sunrealtype** Hes;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform the product :math:`Av`,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``V, q``   - vectors used for workspace by the Arnoldi algorithm.

* ``kry_dim`` - dimension of Krylov subspaces (default is 3),

* ``num_warmups`` - number of preprocessing iterations (default is 100),

* ``num_iters`` - number of iterations (preprocessing and estimation) in the
  last :c:func:`SUNDomEigEstimator_Estimate` call,

* ``num_ATimes`` - number of calls to the ``ATimes`` function,

* ``LAPACK_A, LAPACK_wr, LAPACK_wi, LAPACK_work`` - ``sunrealtype`` used for workspace by LAPACK,

* ``LAPACK_lwork`` - the size of the ``LAPACK_work`` requested by LAPACK,

* ``LAPACK_arr`` - storage for the estimated dominant eigenvalues,

* ``Hes`` - Hessenberg matrix,


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* SUNDIALS packages will call :c:func:`SUNDomEigEstimator_SetATimes` to supply
  the ``ATimes`` function pointer and the related data ``ATData``.

* In :c:func:`SUNDomEigEstimator_Initialize`, the estimator parameters are
  checked for validity and the remaining Arnoldi estimator memory such as LAPACK
  workspace is allocated.

* In :c:func:`SUNDomEigEstimator_Estimate`, the initial nonzero vector
  :math:`q_0` is preprocessed with some fixed number of Power iterations,

  .. math::

     q_1 = \frac{Aq_0}{||Aq_0||} \quad \cdots \quad q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||},

  (see :c:func:`LSRKStepSetNumDomEigEstInitPreprocessIters` and
  :c:func:`LSRKStepSetNumDomEigEstPreprocessIters` for setting the number of
  preprocessing iterations). Then, the Arnoldi iteration is performed to compute
  the estimate.

The SUNDomEigEstimator_Arnoldi module defines implementations of all dominant
eigenvalue estimator operations listed in :numref:`SUNDomEigEst.API`:

* ``SUNDomEigEstimator_SetATimes_Arnoldi``

* ``SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi``

* ``SUNDomEigEstimator_Initialize_Arnoldi``

* ``SUNDomEigEstimator_Estimate_Arnoldi``

* ``SUNDomEigEstimator_GetNumIters_Arnoldi``

* ``SUNDomEigEstimator_GetNumATimesCalls_Arnoldi``

* ``SUNDomEigEstimator_Write_Arnoldi``

* ``SUNDomEigEstimator_Destroy_Arnoldi``
