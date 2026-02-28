.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ UMBC
                  David J. Gardner @ LLNL
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

.. _ARKODE.Usage.UserSupplied:

User-supplied functions
=============================

The user-supplied functions for ARKODE consist of:

* at least one function :ref:`defining the ODE <ARKODE.Usage.ODERHS>`
  (required),

* a function that
  :ref:`provides the error weight vector <ARKODE.Usage.ErrorWeight>` (optional),

* a function that
  :ref:`provides the residual weight vector <ARKODE.Usage.ResidualWeight>` (optional),

* a function that
  :ref:`handles explicit time step stability <ARKODE.Usage.StabilityFn>` (optional),

* a function that
  :ref:`updates the implicit stage prediction <ARKODE.Usage.StagePredictFn>` (optional),

* a function that
  :ref:`defines auxiliary temporal root-finding problem(s) to solve <ARKODE.Usage.RootfindingFn>` (optional),

* one or two functions that
  :ref:`provide Jacobian-related information <ARKODE.Usage.JacobianFn>`
  for the linear solver, if a component is treated implicitly and a
  Newton-based nonlinear iteration is chosen (optional),

* one or two functions that :ref:`define the preconditioner <ARKODE.Usage.PrecSolveFn>`
  for use in any of the Krylov iterative algorithms, if linear systems of
  equations are to be solved using an iterative method (optional),

* if the problem involves a non-identity mass matrix :math:`M\ne I` with ARKStep:

  * one or two functions that
    :ref:`provide mass-matrix-related information <ARKODE.Usage.MassFn>`
    for the linear and mass matrix solvers (required),

  * one or two functions that
    :ref:`define the mass matrix preconditioner <ARKODE.Usage.MassPrecSolveFn>`
    for use if an iterative mass matrix solver is chosen (optional), and

* a function that
  :ref:`handles vector resizing operations <ARKODE.Usage.VecResizeFn>`, if the
  underlying vector structure supports resizing (as opposed to
  deletion/recreation), and if the user plans to call
  :c:func:`ARKodeResize` (optional).

* MRIStep only: functions to be
  :ref:`called before and after each inner integration <ARKODE.Usage.PreInnerFn>`
  to perform any communication or memory transfers of forcing data supplied
  by the outer integrator to the inner integrator, or state data supplied
  by the inner integrator to the outer integrator.

* if relaxation is enabled (optional), a function that
  :ref:`evaluates the conservative or dissipative function <ARKODE.Usage.RelaxFn>`
  :math:`\xi(y(t))` (required) and a function to
  :ref:`evaluate its Jacobian <ARKODE.Usage.RelaxJacFn>`
  :math:`\xi'(y(t))` (required).


.. _ARKODE.Usage.ODERHS:

ODE right-hand side
-----------------------------

The user must supply at least one function of type :c:type:`ARKRhsFn` to
specify the IVP-defininig right-hand side function(s) when creating the
ARKODE time-stepping module:

.. c:type:: int (*ARKRhsFn)(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)

   These functions compute the ODE right-hand side for a given
   value of the independent variable :math:`t` and state vector :math:`y`.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param ydot: the output vector that forms [a portion of] the ODE RHS :math:`f(t,y)`.
   :param user_data: the `user_data` pointer that was passed to
                     :c:func:`ARKodeSetUserData`.

   :return: An *ARKRhsFn* should return 0 if successful, a positive value if a
            recoverable error occurred (in which case ARKODE will attempt to
            correct), or a negative value if it failed unrecoverably (in which
            case the integration is halted and *ARK_RHSFUNC_FAIL* is returned).

   .. note::

      Allocation of memory for `ydot` is handled within ARKODE.

      The vector *ydot* may be uninitialized on input; it is the user's
      responsibility to fill this entire vector with meaningful values.

      A recoverable failure error return from the *ARKRhsFn* is typically
      used to flag a value of the dependent variable :math:`y` that is
      "illegal" in some way (e.g., negative where only a
      non-negative value is physically meaningful).  If such a return is
      made, ARKODE will attempt to recover (possibly repeating the
      nonlinear iteration, or reducing the step size in ARKodeEvolve)
      in order to avoid this recoverable error return.  There are some
      situations in which recovery is not possible even if the right-hand
      side function returns a recoverable error flag.  One is when this
      occurs at the very first call to the *ARKRhsFn* (in which case
      ARKODE returns *ARK_FIRST_RHSFUNC_ERR*).  Another is when a
      recoverable error is reported by *ARKRhsFn* after the time-stepping
      module completes a successful stage, in which case ARKodeEvolve returns
      *ARK_UNREC_RHSFUNC_ERR*).  Finally, when ARKODE is run in fixed-step
      mode, it may halt on a recoverable error flag that would normally have
      resulted in a stepsize reduction.



.. _ARKODE.Usage.ErrorWeight:

Error weight function
--------------------------------------

As an alternative to providing the relative and absolute tolerances,
the user may provide a function of type :c:type:`ARKEwtFn` to compute a
vector *ewt* containing the weights in the WRMS norm
:math:`\|v\|_{WRMS} = \left(\dfrac{1}{n} \displaystyle \sum_{i=1}^n \left(ewt_i\; v_i\right)^2
\right)^{1/2}`.  These weights will be used in place of those defined
in :numref:`ARKODE.Mathematics.Error.Norm`.



.. c:type:: int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void* user_data)

   This function computes the WRMS error weights for the vector
   :math:`y`.

   :param y: the dependent variable vector at which the weight vector is to be computed.
   :param ewt: the output vector containing the error weights.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData` function

   :return: An *ARKEwtFn* function must return 0 if it
            successfully set the error weights, and -1 otherwise.

   .. note::

      Allocation of memory for *ewt* is handled within ARKODE.

      The error weight vector must have all components positive.  It is
      the user's responsibility to perform this test and return -1 if it
      is not satisfied.



.. _ARKODE.Usage.ResidualWeight:

Residual weight function
------------------------

.. warning::

   The functions in this section are specific to time-stepping modules
   that support non-identity mass matrices.


As an alternative to providing the scalar or vector absolute residual
tolerances (when the IVP units differ from the solution units), the
user may provide a function of type :c:type:`ARKRwtFn` to compute a
vector *rwt* containing the weights in the WRMS norm
:math:`\|v\|_{WRMS} = \left(\dfrac{1}{n} \displaystyle \sum_{i=1}^n \left(rwt_i\; v_i\right)^2
\right)^{1/2}`.  These weights will be used in place of those defined
in :numref:`ARKODE.Mathematics.Error.Norm`.



.. c:type:: int (*ARKRwtFn)(N_Vector y, N_Vector rwt, void* user_data)

   This function computes the WRMS residual weights for the vector
   :math:`y`.

   :param y: the dependent variable vector at which the
             weight vector is to be computed.
   :param rwt: the output vector containing the residual weights.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: An *ARKRwtFn* function must return 0 if it
            successfully set the residual weights, and -1 otherwise.

   .. note::

      Allocation of memory for *rwt* is handled within ARKODE.

      The residual weight vector must have all components positive.  It is
      the user's responsibility to perform this test and return -1 if it
      is not satisfied.



.. _ARKODE.Usage.AdaptivityFn:

Time step adaptivity function
--------------------------------------------------------

.. warning::

   The function in this section is only used in now-deprecated functions
   in ARKStep and ERKStep, and will be removed in a future release.


As an alternative to using one of the built-in time step adaptivity
methods for controlling solution error, the user may provide a
function of type :c:type:`ARKAdaptFn` to compute a target step size
:math:`h` for the next integration step.  These steps should be chosen
such that the error estimate for the next time step remains below 1.



.. c:type:: int (*ARKAdaptFn)(N_Vector y, sunrealtype t, sunrealtype h1, sunrealtype h2, sunrealtype h3, sunrealtype e1, sunrealtype e2, sunrealtype e3, int q, int p, sunrealtype* hnew, void* user_data)

   This function implements a time step adaptivity algorithm
   that chooses :math:`h` to satisfy the error tolerances.

   :param y: the current value of the dependent variable vector.
   :param t: the current value of the independent variable.
   :param h1: the current step size, :math:`t_n - t_{n-1}`.
   :param h2: the previous step size, :math:`t_{n-1} - t_{n-2}`.
   :param h3: the step size :math:`t_{n-2}-t_{n-3}`.
   :param e1: the error estimate from the current step, :math:`n`.
   :param e2: the error estimate from the previous step, :math:`n-1`.
   :param e3: the error estimate from the step :math:`n-2`.
   :param q: the global order of accuracy for the method.
   :param p: the global order of accuracy for the embedded method.
   :param hnew: the output value of the next step size.
   :param user_data: a pointer to user data, the same as the
                     *h_data* parameter that was passed to :c:func:`ARKStepSetAdaptivityFn`
                     or :c:func:`ERKStepSetAdaptivityFn`.

   :return: An *ARKAdaptFn* function should return 0 if it
            successfully set the next step size, and a non-zero value otherwise.

   .. deprecated:: 5.7.0

      Use the SUNAdaptController infrastructure instead (see
      :numref:`SUNAdaptController.Description`).



.. _ARKODE.Usage.StabilityFn:

Explicit stability function
---------------------------

.. warning::

   The functions in this section are specific to time-stepping modules
   that support temporal adaptivity.


A user may supply a function to predict the maximum stable step size
for an explicit portion of their IVP.  While the accuracy-based time step
adaptivity algorithms may be sufficient for retaining a stable solution to
the ODE system, these may be inefficient if the explicit right-hand side
function contains moderately stiff terms.  In this scenario, a user may
provide a function of type :c:type:`ARKExpStabFn`
to provide this stability information to ARKODE.  This function
must set the scalar step size satisfying the stability restriction for
the upcoming time step.  This value will subsequently be bounded by
the user-supplied values for the minimum and maximum allowed time
step, and the accuracy-based time step.



.. c:type:: int (*ARKExpStabFn)(N_Vector y, sunrealtype t, sunrealtype* hstab, void* user_data)

   This function predicts the maximum stable step size for the
   explicit portion of the ODE system.

   :param y: the current value of the dependent variable vector.
   :param t: the current value of the independent variable.
   :param hstab: the output value with the absolute value of the
                 maximum stable step size.
   :param user_data: a pointer to user data, the same as the *estab_data*
                     parameter that was passed to :c:func:`ARKodeSetStabilityFn`.

   :return: An *ARKExpStabFn* function should return 0 if it
            successfully set the upcoming stable step size, and a non-zero
            value otherwise.

   .. note::

      If this function is not supplied, or if it returns
      *hstab* :math:`\le 0.0`, then ARKODE will assume that there is no explicit
      stability restriction on the time step size.




.. _ARKODE.Usage.StagePredictFn:

Implicit stage prediction function
----------------------------------

A user may supply a function to update the prediction for each implicit stage solution.
If supplied, this routine will be called *after* any existing ARKODE predictor
algorithm completes, so that the predictor may be modified by the user as desired.
In this scenario, a user may provide a function of type :c:type:`ARKStagePredictFn`
to provide this implicit predictor to ARKODE.  This function takes as input the
already-predicted implicit stage solution and the corresponding "time" for that prediction;
it then updates the prediction vector as desired.  If the user-supplied routine will
construct a full prediction (and thus the ARKODE prediction is irrelevant), it is
recommended that the user *not* call :c:func:`ARKodeSetPredictorMethod`, thereby leaving
the default trivial predictor in place.



.. c:type:: int (*ARKStagePredictFn)(sunrealtype t, N_Vector zpred, void* user_data)

   This function updates the prediction for the implicit stage solution.

   :param t: the current value of the independent variable containing the
             "time" corresponding to the predicted solution.
   :param zpred: the ARKODE-predicted stage solution on input, and the
                 user-modified predicted stage solution on output.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: An *ARKStagePredictFn* function should return 0 if it
            successfully set the upcoming stable step size, and a non-zero
            value otherwise.

   .. note::

      This may be useful if there are bound constraints on the solution,
      and these should be enforced prior to beginning the nonlinear or linear implicit solver
      algorithm.

      This routine is incompatible with the "minimum correction predictor" -- option 5 to the
      routine :c:func:`ARKodeSetPredictorMethod`.  If both are selected, then ARKODE will
      override its built-in implicit predictor routine to instead use option 0 (trivial predictor).


.. _ARKODE.Usage.RootfindingFn:

Rootfinding function
--------------------------------------

If a rootfinding problem is to be solved during integration of the
ODE system, the user must supply a function of type :c:type:`ARKRootFn`.



.. c:type:: int (*ARKRootFn)(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)

   This function implements a vector-valued function
   :math:`g(t,y)` such that roots are sought for the components
   :math:`g_i(t,y)`, :math:`i=0,\ldots,` *nrtfn*-1.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param gout: the output array, of length *nrtfn*, with components :math:`g_i(t,y)`.
   :param user_data: a pointer to user data, the same as the
                     *user_data* parameter that was passed to the ``SetUserData`` function

   :return: An *ARKRootFn* function should return 0 if successful
            or a non-zero value if an error occurred (in which case the
            integration is halted and ARKODE returns *ARK_RTFUNC_FAIL*).

   .. note::

      Allocation of memory for *gout* is handled within ARKODE.



.. _ARKODE.Usage.JacobianFn:

Jacobian construction
---------------------

If a matrix-based linear solver module is used (i.e., a non-NULL ``SUNMatrix``
object was supplied to :c:func:`ARKodeSetLinearSolver`, the user may provide a
function of type :c:type:`ARKLsJacFn` to provide the Jacobian approximation or
:c:type:`ARKLsLinSysFn` to provide an approximation of the linear system
:math:`\mathcal{A}(t,y) = M(t) - \gamma J(t,y)`.


.. c:type:: int (*ARKLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the Jacobian matrix :math:`J(t,y) =
   \dfrac{\partial f^I}{\partial y}(t,y)` (or an approximation to it).

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector, namely
             the predicted value of :math:`y(t)`.
   :param fy: the current value of the vector :math:`f^I(t,y)`.
   :param Jac: the output Jacobian matrix.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.
   :param tmp*: pointers to memory allocated to
                variables of type ``N_Vector`` which can be used by an
                ARKLsJacFn as temporary storage or work space.

   :return: An *ARKLsJacFn* function should return 0 if successful, a positive
            value if a recoverable error occurred (in which case ARKODE will
            attempt to correct, while ARKLS sets *last_flag* to
            *ARKLS_JACFUNC_RECVR*), or a negative value if it failed
            unrecoverably (in which case the integration is halted,
            :c:func:`ARKodeEvolve` returns *ARK_LSETUP_FAIL* and
            ARKLS sets *last_flag* to *ARKLS_JACFUNC_UNRECVR*).

   .. note::

      Information regarding the specific
      ``SUNMatrix`` structure (e.g.~number of rows, upper/lower
      bandwidth, sparsity type) may be obtained through using the
      implementation-specific ``SUNMatrix`` interface functions
      (see :numref:`SUNMatrix` for details).

      When using a linear solver of type ``SUNLINEARSOLVER_DIRECT``, prior
      to calling the user-supplied Jacobian function, the Jacobian
      matrix :math:`J(t,y)` is zeroed out, so only nonzero elements need
      to be loaded into *Jac*.

      With the default Newton nonlinear solver, each
      call to the user's :c:func:`ARKLsJacFn` function is preceded by a call to the
      implicit :c:func:`ARKRhsFn` user function with the same :math:`(t,y)`
      arguments. Thus, the Jacobian function can use any auxiliary data that is
      computed and saved during the evaluation of :math:`f^I(t,y)`.
      In the case of a user-supplied or external nonlinear solver, this is also
      true if the nonlinear system function is evaluated prior to calling the
      linear solver setup function (see :numref:`SUNNonlinSol.API.SUNSuppliedFn` for more
      information).

      If the user's :c:type:`ARKLsJacFn` function uses difference
      quotient approximations, then it may need to access quantities not
      in the argument list, including the current step size, the
      error weights, etc.  To obtain these, the user will need to add a
      pointer to the ``ark_mem`` structure to their ``user_data``, and
      then use the ``ARKSodeGet*`` functions listed in
      :numref:`ARKODE.Usage.OptionalOutputs`. The unit roundoff can be
      accessed as ``SUN_UNIT_ROUNDOFF``, which is defined in the header
      file ``sundials_types.h``.

      **dense** :math:`J(t,y)`:
      A user-supplied dense Jacobian function must load the
      *N* by *N* dense matrix *Jac* with an approximation to the Jacobian
      matrix :math:`J(t,y)` at the point :math:`(t,y)`. Utility routines
      and accessor macros for the SUNMATRIX_DENSE module are documented
      in :numref:`SUNMatrix.Dense`.

      **banded** :math:`J(t,y)`:
      A user-supplied banded Jacobian function must load the band
      matrix *Jac* with the elements of the Jacobian
      :math:`J(t,y)` at the point :math:`(t,y)`. Utility routines
      and accessor macros for the SUNMATRIX_BAND module are
      documented in :numref:`SUNMatrix.Band`.

      **sparse** :math:`J(t,y)`:
      A user-supplied sparse Jacobian function must load the
      compressed-sparse-column (CSC) or compressed-sparse-row (CSR)
      matrix *Jac* with an approximation to the Jacobian matrix
      :math:`J(t,y)` at the point :math:`(t,y)`.  Storage for *Jac*
      already exists on entry to this function, although the user should
      ensure that sufficient space is allocated in *Jac* to hold the
      nonzero values to be set; if the existing space is insufficient the
      user may reallocate the data and index arrays as needed.  Utility
      routines and accessor macros for the SUNMATRIX_SPARSE type are
      documented in :numref:`SUNMatrix.Sparse`.



.. c:type:: int (*ARKLsLinSysFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A, SUNMatrix M, sunbooleantype jok, sunbooleantype *jcur, sunrealtype gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the linear system matrix :math:`\mathcal{A}(t,y) = M(t) - \gamma J(t,y)` (or
   an approximation to it).

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector, namely the
             predicted value of :math:`y(t)`.
   :param fy: the current value of the vector :math:`f^I(t,y)`.
   :param A: the output linear system matrix.
   :param M: the current mass matrix (this input is ``NULL`` if :math:`M = I`).
   :param jok: is an input flag indicating whether the Jacobian-related data
               needs to be updated. The *jok* argument provides for the reuse of
               Jacobian data. When *jok* = ``SUNFALSE``, the Jacobian-related data
               should be recomputed from scratch. When *jok* = ``SUNTRUE`` the Jacobian
               data, if saved from the previous call to this function, can be reused
               (with the current value of *gamma*). A call with *jok* = ``SUNTRUE`` can
               only occur after a call with *jok* = ``SUNFALSE``.
   :param jcur: is a pointer to a flag which should be set to ``SUNTRUE`` if
                Jacobian data was recomputed, or set to ``SUNFALSE`` if Jacobian data
                was not recomputed, but saved data was still reused.
   :param gamma: the scalar :math:`\gamma` appearing in the Newton system matrix
                 :math:`\mathcal{A}=M(t)-\gamma J(t,y)`.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.
   :param tmp*: pointers to memory allocated to variables of
                            type ``N_Vector`` which can be used by an ARKLsLinSysFn as temporary
                            storage or work space.

   :return: An *ARKLsLinSysFn* function should return 0 if successful, a positive value
            if a recoverable error occurred (in which case ARKODE will attempt to
            correct, while ARKLS sets *last_flag* to *ARKLS_JACFUNC_RECVR*), or a
            negative value if it failed unrecoverably (in which case the integration is
            halted, :c:func:`ARKodeEvolve` returns *ARK_LSETUP_FAIL* and ARKLS sets
            *last_flag* to *ARKLS_JACFUNC_UNRECVR*).



.. _ARKODE.Usage.JTimesFn:

Jacobian-vector product
-----------------------

When using a matrix-free linear solver module for the implicit
stage solves (i.e., a NULL-valued SUNMATRIX argument was supplied to
:c:func:`ARKodeSetLinearSolver`, the user may provide a function
of type :c:type:`ARKLsJacTimesVecFn` in the following form, to compute
matrix-vector products :math:`Jv`. If such a function is not supplied,
the default is a difference quotient approximation to these products.


.. c:type:: int (*ARKLsJacTimesVecFn)(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy, void* user_data, N_Vector tmp)

   This function computes the product :math:`Jv` where :math:`J(t,y) \approx
   \dfrac{\partial f^I}{\partial y}(t,y)` (or an approximation to it).

   :param v: the vector to multiply.
   :param Jv: the output vector computed.
   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param fy: the current value of the vector :math:`f^I(t,y)`.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.
   :param tmp: pointer to memory allocated to a variable of type
               ``N_Vector`` which can be used as temporary storage or work space.

   :return: The value to be returned by the Jacobian-vector product
            function should be 0 if successful. Any other return value will
            result in an unrecoverable error of the generic Krylov solver,
            in which case the integration is halted.

   .. note::

      If the user's :c:type:`ARKLsJacTimesVecFn` function
      uses difference quotient approximations, it may need to access
      quantities not in the argument list.  These include the current
      step size, the error weights, etc.  To obtain these, the
      user will need to add a pointer to the ``ark_mem`` structure to
      their ``user_data``, and then use the ``ARKodeGet*`` functions
      listed in :numref:`ARKODE.Usage.OptionalOutputs`. The unit roundoff
      can be accessed as ``SUN_UNIT_ROUNDOFF``, which is defined in the
      header file ``sundials_types.h``.




.. _ARKODE.Usage.JTSetupFn:

Jacobian-vector product setup
-----------------------------

If the user's Jacobian-times-vector routine requires that any Jacobian-related data
be preprocessed or evaluated, then this needs to be done in a
user-supplied function of type :c:type:`ARKLsJacTimesSetupFn`,
defined as follows:


.. c:type:: int (*ARKLsJacTimesSetupFn)(sunrealtype t, N_Vector y, N_Vector fy, void* user_data)

   This function preprocesses and/or evaluates any Jacobian-related
   data needed by the Jacobian-times-vector routine.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param fy: the current value of the vector :math:`f^I(t,y)`.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: The value to be returned by the Jacobian-vector setup
            function should be 0 if successful, positive for a recoverable
            error (in which case the step will be retried), or negative for an
            unrecoverable error (in which case the integration is halted).

   .. note::

      Each call to the Jacobian-vector setup function is
      preceded by a call to the implicit :c:type:`ARKRhsFn` user
      function with the same :math:`(t,y)` arguments.  Thus, the setup
      function can use any auxiliary data that is computed and saved
      during the evaluation of the implicit ODE right-hand side.

      If the user's :c:type:`ARKLsJacTimesSetupFn` function uses
      difference quotient approximations, it may need to access
      quantities not in the argument list.  These include the current
      step size, the error weights, etc.  To obtain these, the
      user will need to add a pointer to the ``ark_mem`` structure to
      their ``user_data``, and then use the ``ARKodeGet*`` functions
      listed in :numref:`ARKODE.Usage.OptionalOutputs`. The unit roundoff
      can be accessed as ``SUN_UNIT_ROUNDOFF``, which is defined in the
      header file ``sundials_types.h``.




.. _ARKODE.Usage.PrecSolveFn:

Preconditioner solve
--------------------

If a user-supplied preconditioner is to be used with a SUNLinSol
solver module, then the user must provide a function of type
:c:type:`ARKLsPrecSolveFn` to solve the linear system :math:`Pz=r`,
where :math:`P` corresponds to either a left or right
preconditioning matrix.  Here :math:`P` should approximate (at least
crudely) the Newton matrix :math:`\mathcal{A}(t,y)=M(t)-\gamma J(t,y)`,
where :math:`M(t)` is the mass matrix and :math:`J(t,y) = \dfrac{\partial f^I}{\partial
y}(t,y)`  If preconditioning is done on both sides, the product of the two
preconditioner matrices should approximate :math:`\mathcal{A}`.



.. c:type:: int (*ARKLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, void* user_data)

   This function solves the preconditioner system :math:`Pz=r`.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param fy: the current value of the vector :math:`f^I(t,y)`.
   :param r: the right-hand side vector of the linear system.
   :param z: the computed output solution vector.
   :param gamma: the scalar :math:`\gamma` appearing in the Newton
                 matrix given by :math:`\mathcal{A}=M(t)-\gamma J(t,y)`.
   :param delta: an input tolerance to be used if an iterative method
                 is employed in the solution.  In that case, the residual vector
                 :math:`Res = r-Pz` of the system should be made to be less than *delta*
                 in the weighted :math:`l_2` norm, i.e. :math:`\left(\displaystyle \sum_{i=1}^n
                 \left(Res_i * ewt_i\right)^2 \right)^{1/2} < \delta`, where :math:`\delta =`
                 `delta`.  To obtain the ``N_Vector`` *ewt*, call
                 :c:func:`ARKodeGetErrWeights`.
   :param lr: an input flag indicating whether the preconditioner
              solve is to use the left preconditioner (*lr* = 1) or the right
              preconditioner (*lr* = 2).
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: The value to be returned by the preconditioner solve
            function is a flag indicating whether it was successful. This value
            should be 0 if successful, positive for a recoverable error (in
            which case the step will be retried), or negative for an
            unrecoverable error (in which case the integration is halted).




.. _ARKODE.Usage.PrecSetupFn:

Preconditioner setup
--------------------

If the user's preconditioner routine above requires that any data be
preprocessed or evaluated, then these actions need to occur within a
user-supplied function of type :c:type:`ARKLsPrecSetupFn`.


.. c:type:: int (*ARKLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)

   This function preprocesses and/or evaluates Jacobian-related
   data needed by the preconditioner.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param fy: the current value of the vector :math:`f^I(t,y)`.
   :param jok: is an input flag indicating whether the Jacobian-related
               data needs to be updated. The *jok* argument provides for the
               reuse of Jacobian data in the preconditioner solve function. When
               *jok* = ``SUNFALSE``, the Jacobian-related data should be recomputed
               from scratch. When *jok* = ``SUNTRUE`` the Jacobian data, if saved from the
               previous call to this function, can be reused (with the current
               value of *gamma*). A call with *jok* = ``SUNTRUE`` can only occur
               after a call with *jok* = ``SUNFALSE``.
   :param jcurPtr: is a pointer to a flag which should be set to
                   ``SUNTRUE`` if Jacobian data was recomputed, or set to ``SUNFALSE`` if
                   Jacobian data was not recomputed, but saved data was still reused.
   :param gamma: the scalar :math:`\gamma` appearing in the Newton
                 matrix given by :math:`\mathcal{A}=M(t)-\gamma J(t,y)`.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: The value to be returned by the preconditioner setup
            function is a flag indicating whether it was successful. This value
            should be 0 if successful, positive for a recoverable error (in
            which case the step will be retried), or negative for an
            unrecoverable error (in which case the integration is halted).

   .. note::

      The operations performed by this function might include
      forming a crude approximate Jacobian, and performing an LU
      factorization of the resulting approximation to :math:`\mathcal{A} = M(t) -
      \gamma J(t,y)`.

      With the default nonlinear solver (the native SUNDIALS Newton method), each
      call to the preconditioner setup function is preceded by a call to the
      implicit :c:type:`ARKRhsFn` user function with the same :math:`(t,y)`
      arguments.  Thus, the preconditioner setup function can use any auxiliary
      data that is computed and saved during the evaluation of the implicit ODE
      right-hand side. In the case of a user-supplied or external nonlinear solver,
      this is also true if the nonlinear system function is evaluated prior to
      calling the linear solver setup function (see
      :numref:`SUNNonlinSol.API.SUNSuppliedFn` for more information).

      This function is not called in advance of every call to the
      preconditioner solve function, but rather is called only as often
      as needed to achieve convergence in the Newton iteration.

      If the user's :c:type:`ARKLsPrecSetupFn` function uses
      difference quotient approximations, it may need to access
      quantities not in the call list. These include the current step
      size, the error weights, etc.  To obtain these, the user will need
      to add a pointer to the ``ark_mem`` structure to their
      ``user_data``, and then use the ``ARKodeGet*`` functions listed in
      :numref:`ARKODE.Usage.OptionalOutputs`. The unit roundoff can be
      accessed as ``SUN_UNIT_ROUNDOFF``, which is defined in the header
      file ``sundials_types.h``.



.. _ARKODE.Usage.MassFn:

Mass matrix construction
------------------------

For problems involving a non-identity mass matrix, if a matrix-based
mass-matrix linear solver is used (i.e., a non-NULL SUNMATRIX was
supplied to :c:func:`ARKodeSetMassLinearSolver`, the
user must provide a function of type :c:type:`ARKLsMassFn` to provide
the mass matrix approximation.



.. c:type:: int (*ARKLsMassFn)(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the mass matrix :math:`M(t)` (or an approximation to it).

   :param t: the current value of the independent variable.
   :param M: the output mass matrix.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.
   :param tmp1*: pointers to memory allocated to
                 variables of type ``N_Vector`` which can be used by an
                 ARKLsMassFn as temporary storage or work space.

   :return: An *ARKLsMassFn* function should return 0 if successful, or a
            negative value if it failed unrecoverably (in which case the
            integration is halted, :c:func:`ARKodeEvolve` returns
            *ARK_MASSSETUP_FAIL* and ARKLS sets *last_flag* to
            *ARKLS_MASSFUNC_UNRECVR*).

   .. note::

      Information regarding the structure of the specific
      ``SUNMatrix`` structure (e.g.~number of rows, upper/lower
      bandwidth, sparsity type) may be obtained through using the
      implementation-specific ``SUNMatrix`` interface functions
      (see :numref:`SUNMatrix` for details).

      Prior to calling the user-supplied mass matrix function, the mass
      matrix :math:`M(t)` is zeroed out, so only nonzero elements need to
      be loaded into *M*.

      **dense** :math:`M(t)`:
      A user-supplied dense mass matrix function must load the *N* by *N*
      dense matrix *M* with an approximation to the mass matrix
      :math:`M(t)`. Utility routines and accessor macros for the
      SUNMATRIX_DENSE module are documented in :numref:`SUNMatrix.Dense`.

      **banded** :math:`M(t)`:
      A user-supplied banded mass matrix function must load the band
      matrix *M* with the elements of the mass matrix :math:`M(t)`.
      Utility routines and accessor macros for the SUNMATRIX_BAND module
      are documented in :numref:`SUNMatrix.Band`.

      **sparse** :math:`M(t)`:
      A user-supplied sparse mass matrix function must load the
      compressed-sparse-column (CSR) or compressed-sparse-row (CSR)
      matrix *M* with an approximation to the mass matrix :math:`M(t)`.
      Storage for *M* already exists on entry to this function, although
      the user should ensure that sufficient space is allocated in *M*
      to hold the nonzero values to be set; if the existing space is
      insufficient the user may reallocate the data and row index arrays
      as needed.  Utility routines and accessor macros for the
      SUNMATRIX_SPARSE type are documented in :numref:`SUNMatrix.Sparse`.




.. _ARKODE.Usage.MTimesFn:

Mass matrix-vector product
--------------------------

For problems involving a non-identity mass matrix, if a matrix-free
linear solver is to be used for mass-matrix linear systems (i.e., a
NULL-valued SUNMATRIX argument was supplied to
:c:func:`ARKodeSetMassLinearSolver` in
:numref:`ARKODE.Usage.Skeleton`), the user *must* provide a
function of type :c:type:`ARKLsMassTimesVecFn` in the following form, to
compute matrix-vector products :math:`M(t)\, v`.



.. c:type:: int (*ARKLsMassTimesVecFn)(N_Vector v, N_Vector Mv, sunrealtype t, void* mtimes_data)

   This function computes the product :math:`M(t)\, v` (or an approximation to it).

   :param v: the vector to multiply.
   :param Mv: the output vector computed.
   :param t: the current value of the independent variable.
   :param mtimes_data: a pointer to user data, the same as the *mtimes_data*
                       parameter that was passed to :c:func:`ARKodeSetMassTimes`.

   :return: The value to be returned by the mass-matrix-vector product
            function should be 0 if successful. Any other return value will
            result in an unrecoverable error of the generic Krylov solver,
            in which case the integration is halted.



.. _ARKODE.Usage.MTSetupFn:

Mass matrix-vector product setup
--------------------------------

For problems involving a non-identity mass matrix and a matrix-free linear
solver, if the user's mass-matrix-times-vector routine requires that any
mass matrix-related data be preprocessed or evaluated, then this needs to
be done in a user-supplied function of type
:c:type:`ARKLsMassTimesSetupFn`, defined as follows:



.. c:type:: int (*ARKLsMassTimesSetupFn)(sunrealtype t, void* mtimes_data)

   This function preprocesses and/or evaluates any mass-matrix-related
   data needed by the mass-matrix-times-vector routine.

   :param t: the current value of the independent variable.
   :param mtimes_data: a pointer to user data, the same as the *mtimes_data*
                       parameter that was passed to :c:func:`ARKodeSetMassTimes`.

   :return: The value to be returned by the mass-matrix-vector setup
            function should be 0 if successful. Any other return value will
            result in an unrecoverable error of the ARKLS mass matrix solver
            interface, in which case the integration is halted.



.. _ARKODE.Usage.MassPrecSolveFn:

Mass matrix preconditioner solve
--------------------------------

For problems involving a non-identity mass matrix and an iterative linear
solver, if a user-supplied preconditioner is to be used with a SUNLINEAR
solver module for mass matrix linear systems, then the user must
provide a function of type :c:type:`ARKLsMassPrecSolveFn` to solve the
linear system :math:`Pz=r`, where :math:`P` may be either a left or right
preconditioning matrix.  Here :math:`P` should approximate (at least
crudely) the mass matrix :math:`M(t)`.  If preconditioning is done on
both sides, the product of the two preconditioner matrices should
approximate :math:`M(t)`.


.. c:type:: int (*ARKLsMassPrecSolveFn)(sunrealtype t, N_Vector r, N_Vector z, sunrealtype delta, int lr, void* user_data)

   This function solves the preconditioner system :math:`Pz=r`.

   :param t: the current value of the independent variable.
   :param r: the right-hand side vector of the linear system.
   :param z: the computed output solution vector.
   :param delta: an input tolerance to be used if an iterative method
                 is employed in the solution.  In that case, the residual vector
                 :math:`Res = r-Pz` of the system should be made to be less than *delta*
                 in the weighted :math:`l_2` norm, i.e. :math:`\left(\displaystyle \sum_{i=1}^n
                 \left(Res_i * ewt_i\right)^2 \right)^{1/2} < \delta`, where :math:`\delta =`
                 *delta*.  To obtain the ``N_Vector`` *ewt*, call
                 :c:func:`ARKodeGetErrWeights`.
   :param lr: an input flag indicating whether the preconditioner
              solve is to use the left preconditioner (*lr* = 1) or the right
              preconditioner (*lr* = 2).
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: The value to be returned by the preconditioner solve
            function is a flag indicating whether it was successful. This value
            should be 0 if successful, positive for a recoverable error (in
            which case the step will be retried), or negative for an
            unrecoverable error (in which case the integration is halted).




.. _ARKODE.Usage.MassPrecSetupFn:

Mass matrix preconditioner setup
--------------------------------

For problems involving a non-identity mass matrix and an iterative linear
solver, if the user's mass matrix preconditioner above requires that any
problem data be preprocessed or evaluated, then these actions need to
occur within a user-supplied function of type
:c:type:`ARKLsMassPrecSetupFn`.



.. c:type:: int (*ARKLsMassPrecSetupFn)(sunrealtype t, void* user_data)

   This function preprocesses and/or evaluates mass-matrix-related
   data needed by the preconditioner.

   :param t: the current value of the independent variable.
   :param user_data: a pointer to user data, the same as the *user_data*
                     parameter that was passed to :c:func:`ARKodeSetUserData`.

   :return: The value to be returned by the mass matrix preconditioner setup
            function is a flag indicating whether it was successful. This value
            should be 0 if successful, positive for a recoverable error (in
            which case the step will be retried), or negative for an
            unrecoverable error (in which case the integration is halted).

   .. note::

      The operations performed by this function might include
      forming a mass matrix and performing an incomplete
      factorization of the result.  Although such operations would
      typically be performed only once at the beginning of a simulation,
      these may be required if the mass matrix can change as a function
      of time.

      If both this function and a :c:type:`ARKLsMassTimesSetupFn` are
      supplied, all calls to this function will be preceded by a call to
      the :c:type:`ARKLsMassTimesSetupFn`, so any setup performed
      there may be reused.


.. _ARKODE.Usage.VecResizeFn:

Vector resize function
--------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatial adaptivity in a
PDE simulation), the ARKODE integrator may be "resized" between
integration steps, through calls to the :c:func:`ARKodeResize`
function. Typically, when performing adaptive simulations the solution
is stored in a customized user-supplied data structure, to enable
adaptivity without repeated allocation/deallocation of memory.  In
these scenarios, it is recommended that the user supply a customized
vector kernel to interface between SUNDIALS and their problem-specific
data structure.  If this vector kernel includes a function of type
:c:type:`ARKVecResizeFn` to resize a given vector implementation, then
this function may be supplied to :c:func:`ARKodeResize` so that all
internal ARKODE vectors may be resized, instead of deleting and
re-creating them at each call.  This resize function should have the
following form:


.. c:type:: int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate, void* user_data)

   This function resizes the vector *y* to match the dimensions of the
   supplied vector, *ytemplate*.

   :param y: the vector to resize.
   :param ytemplate: a vector of the desired size.
   :param user_data: a pointer to user data, the same as the *resize_data*
                     parameter that was passed to :c:func:`ARKodeResize`.

   :return: An *ARKVecResizeFn* function should return 0 if it successfully
            resizes the vector *y*, and a non-zero value otherwise.

   .. note::

      If this function is not supplied, then ARKODE will
      instead destroy the vector *y* and clone a new vector *y* off of
      *ytemplate*.




.. _ARKODE.Usage.PreInnerFn:

Pre inner integrator communication function (MRIStep only)
----------------------------------------------------------

The user may supply a function of type :c:type:`MRIStepPreInnerFn` that will be
called *before* each inner integration to perform any communication or
memory transfers of forcing data supplied by the outer integrator to the inner
integrator for the inner integration.


.. c:type:: int (*MRIStepPreInnerFn)(sunrealtype t, N_Vector* f, int num_vecs, void* user_data)

   :param t: the current value of the independent variable.
   :param f: an ``N_Vector`` array of outer forcing vectors.
   :param num_vecs: the number of vectors in the ``N_Vector`` array.
   :param user_data: the `user_data` pointer that was passed to
                     :c:func:`MRIStepSetUserData`.

   :return: An *MRIStepPreInnerFn* function should return 0 if successful, a positive value
            if a recoverable error occurred, or a negative value if an unrecoverable
            error occurred. As the MRIStep module only supports fixed step sizes at this
            time any non-zero return value will halt the integration.

   .. note::

      In a heterogeneous computing environment if any data copies between the host
      and device vector data are necessary, this is where that should occur.


.. _ARKODE.Usage.PostInnerFn:

Post inner integrator communication function (MRIStep only)
-----------------------------------------------------------

The user may supply a function of type :c:type:`MRIStepPostInnerFn` that will be
called *after* each inner integration to perform any communication or
memory transfers of state data supplied by the inner integrator to the
outer integrator for the outer integration.


.. c:type:: int (*MRIStepPostInnerFn)(sunrealtype t, N_Vector y, void* user_data)

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param user_data: the ``user_data`` pointer that was passed to
                     :c:func:`MRIStepSetUserData`.

   :return:  An :c:func:`MRIStepPostInnerFn` function should return 0 if successful, a
             positive value if a recoverable error occurred, or a negative value if an
             unrecoverable error occurred. As the MRIStep module only supports fixed step
             sizes at this time any non-zero return value will halt the integration.

   .. note::

      In a heterogeneous computing environment if any data copies between the host
      and device vector data are necessary, this is where that should occur.


.. _ARKODE.Usage.RelaxFn:

Relaxation function
-------------------

.. c:type:: int (*ARKRelaxFn)(N_Vector y, sunrealtype* r, void* user_data)

   When applying relaxation, an :c:func:`ARKRelaxFn` function is required to
   compute the conservative or dissipative function :math:`\xi(y)`.

   :param y: the current value of the dependent variable vector.
   :param r: the value of :math:`\xi(y)`.
   :param user_data: the ``user_data`` pointer that was passed to
                     :c:func:`ARKodeSetUserData`.

   :return: An :c:func:`ARKRelaxFn` function should return 0 if successful, a positive
            value if a recoverable error occurred, or a negative value if an
            unrecoverable error occurred. If a recoverable error occurs, the step size
            will be reduced and the step repeated.

.. _ARKODE.Usage.RelaxJacFn:

Relaxation Jacobian function
----------------------------

.. c:type:: int (*ARKRelaxJacFn)(N_Vector y, N_Vector J, void* user_data);

   When applying relaxation, an :c:func:`ARKRelaxJacFn` function is required to
   compute the Jacobian :math:`\xi'(y)` of the :c:func:`ARKRelaxFn`
   :math:`\xi(y)`.

   :param y: the current value of the dependent variable vector.
   :param J: the Jacobian vector :math:`\xi'(y)`.
   :param user_data: the ``user_data`` pointer that was passed to
                     :c:func:`ARKodeSetUserData`.

   :return: An :c:func:`ARKRelaxJacFn` function should return 0 if successful, a
            positive value if a recoverable error occurred, or a negative value if an
            unrecoverable error occurred. If a recoverable error occurs, the step size
            will be reduced and the step repeated.
