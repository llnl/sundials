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

.. _SUNDomEigEst.Power:

The SUNDomEigEstimator_Power Module
======================================

.. versionadded:: 7.5.0

The SUNDomEigEstimator_Power implementation of the :c:type:`SUNDomEigEstimator`
class performs the Power Iteration (PI) method :cite:p:`vonmises29`; this is an
iterative dominant eigenvalue estimator that is designed to be compatible with
any ``N_Vector`` implementation that supports a minimal subset of operations
(:c:func:`N_VClone()`, :c:func:`N_VDotProd()`, :c:func:`N_VScale()`, and
:c:func:`N_VDestroy()`).

Power iteration is useful for large, sparse matrices whose dominant eigenvalue
is real-valued and has algebraic multiplicity one. The algorithm starts with a non-zero
vector :math:`\mathbf{v}_{0}`.  It then  iteratively updates this via

.. math::

    \mathbf{v}_{k+1} = \frac{A \mathbf{v}_k}{\|A \mathbf{v}_k\|},

where :math:`\| \cdot \|` denotes the Euclidean norm.  Over successive iterations,
:math:`\mathbf{v}_k` converges to the eigenvector corresponding to
the dominant eigenvalue of :math:`A`.  At each step, the corresponding eigenvalue
can be approximated using the Rayleigh quotient

.. math::

    \lambda_k = \frac{\mathbf{v}_k^T A \mathbf{v}_k}{\|\mathbf{v}_k\|^2}.

The iteration continues until the two successive eigenvalue approximations are
relatively close enough to one another.  That is, for some :ref:`relative tolerance <pi_rel_tol>`.

Power iteration works for the matrices that have a **real** dominant eigenvalue.
If it is strictly greater than all others (in magnitude), convergence is guaranteed.
The speed of convergence depends on the ratios of the magnitude of the first two dominant eigenvalues.

The matrix :math:`A` is not required explicitly; only a routine that provides
the matrix-vector product :math:`Av` is required.  Also, PI requires a fixed
amount of memory regardless of the number of iterations.


.. _SUNDomEigEst.Power.Usage:

SUNDomEigEstimator_Power Usage
------------------------------

To use SUNDomEigEstimator_Arnoldi include the header file
``sundomeigest/sundomeigest_power.h``, and link to the library
``libsundials_sundomeigestpower``.

The module SUNDomEigEstimator_Power provides the following user-callable
routines:


.. c:function:: SUNDomEigEstimator SUNDomEigEstimator_Power(N_Vector q, long int max_iters, sunrealtype rel_tol, SUNContext sunctx)

   This constructor function creates and allocates memory for the Power
   iteration implementation of a :c:type:`SUNDomEigEstimator`.

   Consistency checks are performed to ensure the input vector is non-zero and
   supplies the necessary operations.

   :param q: the initial guess for the dominant eigenvector; this should not
             be a non-dominant eigenvector of the Jacobian.
   :param max_iters: maximum number of iterations (default 100). Supplying a
                     value :math:`\leq 0` will result in using the default
                     value. Although this default number is not high for large
                     matrices, it is reasonable since (1) most solvers do not
                     need too tight tolerances and consider a safety factor,
                     and (2) an early (less costly) termination will be a good
                     indicator whether the power iteration is compatible.
   :param rel_tol: relative tolerance for convergence checks (default 0.005). A
                   value :math:`\leq 0` will result in the default value. The
                   default has been found to small enough for many internal
                   applications.
   :param sunctx: the :c:type:`SUNContext` object.

   :returns: If successful, a :c:type:`SUNDomEigEstimator` otherwise ``NULL``.

   .. note::

      When used in a time-dependent context, the initial guess supplied to the
      constructor, ``q``, is used only for the first
      :c:func:`SUNDomEigEstimator_Estimate` call and is overwritten with the
      result of the next to last Power iteration from the most recent
      :c:func:`SUNDomEigEstimator_Estimate` call. This new value is used as the
      initial guess for subsequent estimates.

      The initial guess can be reset with
      :c:func:`SUNDomEigEstimator_SetInitialGuess`.


.. _SUNDomEigEst.Power.Description:

SUNDomEigEstimator_Power Description
------------------------------------


The SUNDomEigEstimator_Power module defines the *content* field of a
``SUNDomEigEstimator`` to be the following structure:

.. code-block:: c

   struct SUNDomEigEstimatorContent_Power_ {
     SUNATimesFn ATimes;
     void* ATdata;
     N_Vector* V;
     N_Vector q;
     int num_warmups;
     long int max_iters;
     long int num_iters;
     long int num_ATimes;
     sunrealtype rel_tol;
     sunrealtype res;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform the product :math:`Av`,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``V, q``   - ``N_Vector`` used for workspace by the PI algorithm.

* ``num_warmups`` - number of preprocessing iterations (default is 100),

* ``max_iters`` - maximum number of iterations (default is 100),

* ``num_iters`` - number of iterations (preprocessing and estimation) in the
  last :c:func:`SUNDomEigEstimator_Estimate` call,

* ``num_ATimes`` - number of calls to the ``ATimes`` function,

* ``rel_tol`` - relative tolerance for the convergence criteria (default is 0.005),

* ``res`` - the residual from the last :c:func:`SUNDomEigEstimator_Estimate`
  call.


This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* SUNDIALS packages will call :c:func:`SUNDomEigEstimator_SetATimes` to supply
  the ``ATimes`` function pointer and the related data ``ATData``.

* In :c:func:`SUNDomEigEstimator_Initialize`, the estimator parameters are
  checked for validity and the initial eigenvector is normalized.

* In :c:func:`SUNDomEigEstimator_Estimate`, the initial nonzero vector
  :math:`q_0` is preprocessed with some fixed number of Power iterations,

  .. math::

     q_1 = \frac{Aq_0}{||Aq_0||} \quad \cdots \quad q_k = \frac{Aq_{k-1}}{||Aq_{k-1}||},

  (see :c:func:`LSRKStepSetNumDomEigEstInitPreprocessIters` and
  :c:func:`LSRKStepSetNumDomEigEstPreprocessIters` for setting the number of
  preprocessing iterations) before computing the estimate.

The SUNDomEigEstimator_Power module defines implementations of all dominant
eigenvalue estimator operations listed in :numref:`SUNDomEigEst.API`:

* ``SUNDomEigEstimator_SetATimes_Power``

* ``SUNDomEigEstimator_SetMaxIters_Power``

* ``SUNDomEigEstimator_SetNumPreprocessIters_Power``

* ``SUNDomEigEstimator_SetRelTol_Power``

* ``SUNDomEigEstimator_Initialize_Power``

* ``SUNDomEigEstimator_Estimate_Power``

* ``SUNDomEigEstimator_SetInitialGuess_Power``

* ``SUNDomEigEstimator_GetRes_Power``

* ``SUNDomEigEstimator_GetNumIters_Power``

* ``SUNDomEigEstimator_GetNumATimesCalls_Power``

* ``SUNDomEigEstimator_Write_Power``

* ``SUNDomEigEstimator_Destroy_Power``
