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
has algebraic multiplicity one, or if the dominant eigenvalues are a complex 
conjugate pair, then that pair has algebraic multiplicity one.  The algorithm
starts with a non-zero vector :math:`\mathbf{v}_{0}`.  It then  iteratively 
updates this via

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

The default Power iteration implementation estimates complex-valued dominant
eigenvalues with real arithmetic.  After the iteration phase, a postprocessing 
step is performed using the two most recent iterate vectors (approximations of 
the dominant eigenvector).  These vectors are used to construct a :math:`2 \times 2` projection
of the original matrix.

If the two iterates are (numerically) linearly dependent, this indicates
convergence to a one-dimensional invariant subspace, consistent with a
real-valued dominant eigenvalue.  In this case, the dominant eigenvalue
estimate is taken as the Rayleigh quotient of the final iterate.

If the iterates are not linearly dependent, they span a two-dimensional
subspace.  A :math:`2 \times 2` projection of the original matrix onto this subspace is
constructed, and the eigenvalues of this projected matrix are used as the
dominant eigenvalue estimates.  This allows the method to capture complex
conjugate dominant eigenvalue pairs.

An option is also provided to estimate only a real-valued dominant
eigenvalue.  In this mode, the :math:`2 \times 2` projection step is skipped and the
Rayleigh quotient of the final iterate is returned directly.

If the dominant eigenvalue is strictly greater than all others (in magnitude), 
convergence is guaranteed.  The speed of convergence depends on the ratios of 
the magnitude of the first two dominant eigenvalues.

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


.. c:function:: SUNErrCode SUNDomEigEstimator_SetIsReal_Power(SUNDomEigEstimator DEE, sunbooleantype real)

   This routine informs the Power iteration that the dominant eigenvalue is 
   real-valued, so that the complex projection described in Section 
   :numref:`SUNDomEigEst.Power` can be omitted.

   :param DEE: the dominant eigenvalue estimator object.
   :param real: flag indicating that the dominant eigenvalue is real-valued.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an appropriate error code. 

   .. note::

      No matter the value of ``real``, the convergence criterion is based on the relative change in the 
      magnitude of successive eigenvalue estimates (with tolerance set using 
      :c:func:`SUNDomEigEstimator_SetRelTol`).  If ``real`` is ``SUNTRUE``, then the final Power 
      iteration estimate is returned.  Otherwise, a postprocessing step is performed to compute 
      the complex-valued dominant eigenvalue estimate. 

      The default value is ``SUNFALSE``.


.. c:function:: SUNErrCode SUNDomEigEstimator_SetRelTol_Power(SUNDomEigEstimator DEE, sunrealtype rel_tol)

   This routine sets the relative tolerance used by the Power iteration.

   :param DEE: the dominant eigenvalue estimator object.
   :param rel_tol: requested relative tolerance.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an appropriate error code.

   .. note::

      In the Power implementation, ``rel_tol`` is used in two ways.

      First, it is used in the convergence test for successive dominant
      eigenvalue estimates:

      .. math::

         \left|\lambda_{k} - \lambda_{k-1}\right|
         \le \mathtt{rel\_tol} \cdot |\lambda_{k}|.

      Second, when estimating complex-valued dominant eigenvalues with real
      arithmetic, it is used to determine whether the two most recent iterates
      are numerically linearly dependent. In this case,

      .. math::

         \mathtt{gram\_det\_tol} = 10 \cdot \max\left(\varepsilon,\; \mathtt{rel\_tol}\right),

      where :math:`\varepsilon` denotes machine precision. If the determinant
      of the :math:`2 \times 2` Gram matrix formed from the two most recent iterates is less
      than or equal to ``gram_det_tol``, then the dominant eigenvalue is treated
      as real.

      To avoid masking a small imaginary part of the dominant eigenvalue,
      ``rel_tol`` should not be chosen too large. In practice, the smallest
      reliably detectable imaginary part satisfies

      .. math::

         |\beta| \gtrsim \mathcal{O}(\mathtt{rel\_tol}),

      so to resolve an expected imaginary part of magnitude :math:`|\beta|`, it
      is recommended to choose

      .. math::

         \mathtt{rel\_tol} \ll |\beta|.

      Choosing a smaller relative tolerance improves the ability to detect
      weakly complex eigenvalues, but may increase computational cost.

      Acceptable inputs satisfy :math:`0 < \mathtt{rel\_tol} < 1`.
      Values outside this range reset to the default value ``0.005``.


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
     N_Vector Av;
     N_Vector v_prev;
     N_Vector rhs_linY;
     N_Vector Fy;
     N_Vector work;
     int num_warmups;
     long int max_iters;
     long int num_iters;
     long int num_ATimes;
     sunrealtype rhs_linT;
     sunrealtype rel_tol;
     sunrealtype res;
     SUNRhsFn rhsfn;
     void* rhs_data;
     long int nfevals;
     sunbooleantype is_complex;
   };


These entries of the *content* field contain the following
information:

* ``ATimes`` - function pointer to perform the product :math:`Av`,

* ``ATdata`` - pointer to structure for ``ATimes``,

* ``V, Av, v_prev, Fy, work``   - ``N_Vector`` used for workspace.

* ``num_warmups`` - number of preprocessing iterations (default is 100),

* ``max_iters`` - maximum number of iterations (default is 100),

* ``num_iters`` - number of iterations (preprocessing and estimation) in the
  last :c:func:`SUNDomEigEstimator_Estimate` call,

* ``num_ATimes`` - number of calls to the ``ATimes`` function,

* ``rhs_linY`` - state vector for linearization point,

* ``rhs_linT`` - time value for linearization point,

* ``rel_tol`` - relative tolerance for the convergence criteria (default is 0.005),

* ``res`` - the residual from the last :c:func:`SUNDomEigEstimator_Estimate`
  call.

* ``rhsfn`` - user provided RHS function,

* ``rhs_data`` - pointer to the data structure for ``rhsfn``,

* ``nfevals`` - number of RHS evaluations,

* ``is_complex`` - flag indicating whether the dominant eigenvalue is 
  complex-valued (default is ``SUNTRUE``).

This estimator is constructed to perform the following operations:

* During construction all ``N_Vector`` estimator data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default generic estimator parameters are set.

* User-facing "set" routines may be called to modify default
  estimator parameters.

* SUNDIALS packages will call :c:func:`SUNDomEigEstimator_SetATimes` to supply
  the ``ATimes`` function pointer and the related data ``ATData``. Or, the user
  may call :c:func:`SUNDomEigEstimator_SetRhs` to supply the RHS function and
  related data. This approach internally constructs an ``ATimes`` function that
  uses the RHS function to compute the matrix-vector product :math:`Av` for
  the Jacobian of the RHS function.

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

* ``SUNDomEigEstimator_SetRhs_Power``

* ``SUNDomEigEstimator_SetRhsLinearizationPoint_Power``

* ``SUNDomEigEstimator_SetNumPreprocessIters_Power``

* ``SUNDomEigEstimator_SetMaxIters_Power``

* ``SUNDomEigEstimator_SetRelTol_Power``

* ``SUNDomEigEstimator_SetInitialGuess_Power``

* ``SUNDomEigEstimator_SetIsReal_Power``

* ``SUNDomEigEstimator_Initialize_Power``

* ``SUNDomEigEstimator_Estimate_Power``

* ``SUNDomEigEstimator_GetRes_Power``

* ``SUNDomEigEstimator_GetNumIters_Power``

* ``SUNDomEigEstimator_GetNumRhsEvals_Power``

* ``SUNDomEigEstimator_GetNumATimesCalls_Power``

* ``SUNDomEigEstimator_Write_Power``

* ``SUNDomEigEstimator_Destroy_Power``
