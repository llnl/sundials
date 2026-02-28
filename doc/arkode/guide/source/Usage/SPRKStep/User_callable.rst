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

.. _ARKODE.Usage.SPRKStep.UserCallable:

SPRKStep User-callable functions
==================================

This section describes the SPRKStep-specific functions that may be called
by the user to setup and then solve an IVP using the SPRKStep time-stepping
module. The large majority of these routines merely wrap :ref:`underlying
ARKODE functions <ARKODE.Usage.UserCallable>`, and are now deprecated
-- each of these are clearly marked.  However, some
of these user-callable functions are specific to SPRKStep, as explained
below.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
SPRKStep supports only the basic set of user-callable functions, and
does not support any of the restricted groups (time adaptivity, implicit
solvers, etc.).

SPRKStep does not have forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and :c:func:`ARKStepCreateMRIStepInnerStepper`
for additional details.


.. _ARKODE.Usage.SPRKStep.Initialization:

SPRKStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,\
                                     N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the SPRKStep time-stepping module in ARKODE.

   :param f1: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f_1(t,q) = -\frac{\partial V(t,q)}{\partial q}`
   :param f2: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f_2(t,p) = \frac{\partial T(t,p)}{\partial p}`
   :param t0: the initial value of :math:`t`
   :param y0: the initial condition vector :math:`y(t_0)`
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :returns: If successful, a pointer to initialized problem memory of type
             ``void*``, to be passed to all user-facing SPRKStep routines listed
             below.  If unsuccessful, a ``NULL`` pointer will be returned, and
             an error message will be printed to ``stderr``.

   .. warning::

      SPRKStep requires a partitioned problem where ``f1`` should only modify
      the q variables and ``f2`` should only modify the p variables (or vice
      versa). However, the vector passed to these functions is the full vector
      with both p and q. The ordering of the variables is determined implicitly
      by the user when they set the initial conditions.


.. c:function:: void SPRKStepFree(void** arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`SPRKStepCreate`.

   :param arkode_mem: pointer to the SPRKStep memory block.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeFree` instead.


.. _ARKODE.Usage.SPRKStep.RootFinding:

Rootfinding initialization function
--------------------------------------

.. c:function:: int SPRKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system. It *must* be called after
   :c:func:`SPRKStepCreate`, and before :c:func:`SPRKStepEvolve()`.

   To disable the rootfinding feature after it has already been initialized, or
   to free memory associated with SPRKStep's rootfinding module, call
   :c:func:`SPRKStepRootInit` with *nrtfn = 0*.

   Similarly, if a new IVP is to be solved with a call to
   :c:func:`SPRKStepReInit()`, where the new IVP has no rootfinding problem but
   the prior one did, then call :c:func:`SPRKStepRootInit` with *nrtfn = 0*.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nrtfn: number of functions :math:`g_i`, an integer :math:`\ge` 0.
   :param g: name of user-supplied function, of type :c:func:`ARKRootFn()`,
      defining the functions :math:`g_i` whose roots are sought.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL: if there was a memory allocation failure
   :retval ARK_ILL_INPUT: if *nrtfn* is greater than zero but *g* = ``NULL``.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeRootInit` instead.


.. _ARKODE.Usage.SPRKStep.Integration:

SPRKStep solver function
-------------------------


.. c:function:: int SPRKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout, sunrealtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tout: the next time at which a computed solution is desired.
   :param yout: the computed solution vector.
   :param tret: the time corresponding to *yout* (output).
   :param itask: a flag indicating the job of the solver for the next user step.

                 The *ARK_NORMAL* option causes the solver to take internal
                 steps until it has just overtaken a user-specified output
                 time, *tout*, in the direction of integration,
                 i.e. :math:`t_{n-1} <` *tout* :math:`\le t_{n}` for forward
                 integration, or :math:`t_{n} \le` *tout* :math:`< t_{n-1}` for
                 backward integration.  It will then compute an approximation
                 to the solution :math:`y(tout)` by interpolation (using one
                 of the dense output routines described in
                 :numref:`ARKODE.Mathematics.Interpolation`).

                 The *ARK_ONE_STEP* option tells the solver to only take a
                 single internal step, :math:`y_{n-1} \to y_{n}`, and return the
                 solution at that point, :math:`y_{n}`, in the vector *yout*.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_ROOT_RETURN: if :c:func:`SPRKStepEvolve()` succeeded, and
                            found one or more roots.  If the number of root functions,
                            *nrtfn*, is greater than 1, call
                            :c:func:`SPRKStepGetRootInfo()` to see which :math:`g_i` were
                            found to have a root at (*\*tret*).
   :retval ARK_TSTOP_RETURN: if :c:func:`SPRKStepEvolve()` succeeded and
                             returned at *tstop*.
   :retval ARK_MEM_NULL: if the *arkode_mem* argument was ``NULL``.
   :retval ARK_NO_MALLOC: if *arkode_mem* was not allocated.
   :retval ARK_ILL_INPUT: if one of the inputs to
                          :c:func:`SPRKStepEvolve()` is illegal, or some other
                          input to the solver was either illegal or missing.
                          Details will be provided in the error message. Typical
                          causes of this failure are a root of one of the root
                          functions was found both at a point :math:`t` and also
                          very near :math:`t`.
   :retval ARK_TOO_MUCH_WORK: if the solver took *mxstep* internal steps but
                              could not reach *tout*. The default value for
                              *mxstep* is *MXSTEP_DEFAULT = 500*.
   :retval ARK_ERR_FAILURE: if error test failures occurred either too many
                            times (*ark_maxnef*) during one internal time step
                            or occurred with :math:`|h| = h_{min}`.
   :retval ARK_VECTOROP_ERR: a vector operation error occurred.

   .. note::

      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      :c:func:`SPRKStepCreate`.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable.

      All failure return values are negative and so testing the
      return argument for negative values will trap all
      :c:func:`SPRKStepEvolve()` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`SPRKStepSetStopTime()` before the call to
      :c:func:`SPRKStepEvolve()` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`SPRKStepEvolve()`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`SPRKStepSetStopTime()`). Interpolated outputs may or may not
      conserve the Hamiltonian. Our testing has shown that Lagrange
      interpolation typically performs well in this regard, while Hermite
      interpolation does not. As such, SPRKStep uses the Lagrange interpolation
      module by default.

      On any error return in which one or more internal steps were taken
      by :c:func:`SPRKStepEvolve()`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeEvolve` instead.




.. _ARKODE.Usage.SPRKStep.OptionalInputs:

Optional input functions
-------------------------


.. _ARKODE.Usage.SPRKStep.SPRKStepInput:

Optional inputs for SPRKStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int SPRKStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to SPRKStep's original
   default values.

   :param arkode_mem: pointer to the SPRKStep memory block.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. note::

      Does not change problem-defining function pointer *f*
      or the *user_data* pointer.

      Also leaves alone any data structures or options related to
      root-finding (those can be reset using :c:func:`SPRKStepRootInit()`).

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetDefaults` instead.


.. c:function:: int SPRKStepSetInterpolantType(void* arkode_mem, int itype)

   .. deprecated:: 6.1.0

      This function is now a wrapper to :c:func:`ARKodeSetInterpolantType`, see
      the documentation for that function instead.


.. c:function:: int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree)

   Specifies the degree of the polynomial interpolant
   used for dense output (i.e. interpolation of solution output values).
   Allowed values are between 0 and 5.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param degree: requested polynomial degree.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory or interpolation module are ``NULL``
   :retval ARK_INTERP_FAIL: if this is called after :c:func:`SPRKStepEvolve()`
   :retval ARK_ILL_INPUT: if an argument had an illegal value or the
                          interpolation module has already been initialized

   .. note::

      This routine should be called *after* :c:func:`SPRKStepCreate` and *before*
      :c:func:`SPRKStepEvolve()`. After the first call to :c:func:`SPRKStepEvolve()`
      the interpolation degree may not be changed without first calling
      :c:func:`SPRKStepReInit()`.

      If a user calls both this routine and :c:func:`SPRKStepSetInterpolantType()`, then
      :c:func:`SPRKStepSetInterpolantType()` must be called first.

      Since the accuracy of any polynomial interpolant is limited by the
      accuracy of the time-step solutions on which it is based, the *actual*
      polynomial degree that is used by SPRKStep will be the minimum of
      :math:`q-1` and the input *degree*, for :math:`q > 1` where :math:`q` is
      the order of accuracy for the time integration method.

      When `q = 1`, a linear interpolant is the default to ensure values
      obtained by the integrator are returned at the ends of the time interval.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetInterpolantDegree` instead.


.. c:function:: int SPRKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)

   Sets the time step size used within SPRKStep.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param hfixed: value of the fixed step size to use.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetFixedStep` instead.


.. c:function:: int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before SPRKStep
   will return with an error.

   Passing *mxsteps* = 0 results in SPRKStep using the
   default value (500).

   Passing *mxsteps* < 0 disables the test (not recommended).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param mxsteps: maximum allowed number of internal steps.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxNumSteps` instead.


.. c:function:: int SPRKStepSetStopTime(void* arkode_mem, sunrealtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.

   The default is that no stop time is imposed.

   Once the integrator returns at a stop time, any future testing for
   ``tstop`` is disabled (and can be re-enabled only though a new call to
   :c:func:`SPRKStepSetStopTime`).

   A stop time not reached before a call to :c:func:`SPRKStepReInit` or
   :c:func:`SPRKStepReset` will remain active but can be disabled by calling
   :c:func:`SPRKStepClearStopTime`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tstop: stopping time for the integrator.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetStopTime` instead.


.. c:function:: int SPRKStepClearStopTime(void* arkode_mem)

   Disables the stop time set with :c:func:`SPRKStepSetStopTime`.

   The stop time can be re-enabled though a new call to
   :c:func:`SPRKStepSetStopTime`.

   :param arkode_mem: pointer to the SPRKStep memory block.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeClearStopTime` instead.


.. c:function:: int SPRKStepSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* and
   attaches it to the main SPRKStep memory block.

   If specified, the pointer to *user_data* is passed to all
   user-supplied functions for which it is an argument; otherwise
   ``NULL`` is passed.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param user_data: pointer to the user data.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetUserData` instead.


.. _ARKODE.Usage.SPRKStep.SPRKStepMethodInput:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.SPRKStep.SPRKStepMethodInputTable:
.. table:: Optional inputs for IVP method selection

   +-----------------------------+-------------------------------------------+-------------------------------------+
   |       Optional input        |               Function name               |               Default               |
   +=============================+===========================================+=====================================+
   | Set integrator method order | :c:func:`SPRKStepSetOrder()`              | 4                                   |
   +-----------------------------+-------------------------------------------+-------------------------------------+
   | Set SPRK method             | :c:func:`SPRKStepSetMethod()`             | ``ARKODE_SPRK_MCLACHLAN_4_4``       |
   +-----------------------------+-------------------------------------------+-------------------------------------+
   | Set SPRK method by name     | :c:func:`SPRKStepSetMethodName()`         | "ARKODE_SPRK_MCLACHLAN_4_4"         |
   +-----------------------------+-------------------------------------------+-------------------------------------+
   | Use compensated summation   | :c:func:`SPRKStepSetUseCompensatedSums()` | false                               |
   +-----------------------------+-------------------------------------------+-------------------------------------+


.. c:function:: int SPRKStepSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the SPRK integration method.

   The allowed values are :math:`1,2,3,4,5,6,8,10`. Any illegal input will
   result in the default value of 4.

   Since *ord* affects the memory requirements for the internal
   SPRKStep memory block, it cannot be changed after the first call to
   :c:func:`SPRKStepEvolve()`, unless :c:func:`SPRKStepReInit()` is called.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ord: requested order of accuracy.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. warning::

      This overrides any previously set method so it should not be used with
      :c:func:`SPRKStepSetMethod` or :c:func:`SPRKStepSetMethodName`.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetOrder` instead.


.. c:function:: int SPRKStepSetMethod(void* arkode_mem, ARKodeSPRKTable sprk_table)

   Specifies the SPRK method.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param sprk_table: the SPRK method table.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. note::

      No error checking is performed on the coefficients contained in the
      table to ensure its declared order of accuracy.

   .. warning::

      This should not be used with :c:func:`ARKodeSetOrder`.


.. c:function:: int SPRKStepSetMethodName(void* arkode_mem, const char* method)

   Specifies the SPRK method by its name.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param method: the SPRK method name.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.method_name".

   .. warning::

      This should not be used with :c:func:`ARKodeSetOrder`.


.. c:function:: int SPRKStepSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)

   Specifies if :ref:`compensated summation (and the incremental form) <ARKODE.Mathematics.SPRKStep.Compensated>`
   should be used where applicable.

   This increases the computational cost by 2 extra vector operations per stage
   and an additional 5 per time step. It also requires one extra vector to be
   stored. However, it is significantly more robust to roundoff error
   accumulation.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param onoff: should compensated summation be used (1) or not (0)

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. deprecated:: 6.4.0

      Use :c:func:`ARKodeSetUseCompensatedSums` instead.


.. _ARKODE.Usage.SPRKStep.SPRKStepRootfindingInput:


Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int SPRKStepSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.

   The default behavior is to monitor for both zero-crossing directions.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param rootdir: state array of length *nrtfn*, the number of root
      functions :math:`g_i`  (the value of *nrtfn* was supplied in
      the call to :c:func:`SPRKStepRootInit()`).  If ``rootdir[i] ==
      0`` then crossing in either direction for :math:`g_i` should be
      reported.  A value of +1 or -1 indicates that the solver
      should report only zero-crossings where :math:`g_i` is
      increasing or decreasing, respectively.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRootDirection` instead.


.. c:function:: int SPRKStepSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.

   SPRKStep will not report the initial conditions as a possible zero-crossing
   (assuming that one or more components :math:`g_i` are zero at the initial
   time). However, if it appears that some :math:`g_i` is identically zero at
   the initial time (i.e., :math:`g_i` is zero at the initial time *and* after
   the first step), SPRKStep will issue a warning which can be disabled with
   this optional input function.

   :param arkode_mem: pointer to the SPRKStep memory block.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetNoInactiveRootWarn` instead.


.. _ARKODE.Usage.SPRKStep.InterpolatedOutput:

Interpolated output function
--------------------------------

.. c:function:: int SPRKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function :math:`y` at the time *t*,
   i.e., :math:`y^{(k)}(t)`, for values of the independent variable satisfying
   :math:`t_n-h_n \le t \le t_n`, with :math:`t_n` as current internal time
   reached, and :math:`h_n` is the last internal step size successfully used by
   the solver. A user may access the values :math:`t_n` and :math:`h_n` via the
   functions :c:func:`SPRKStepGetCurrentTime()` and
   :c:func:`SPRKStepGetLastStep()`, respectively.

   This routine uses an interpolating polynomial of degree *min(degree, 5)*,
   where *degree* is the argument provided to
   :c:func:`SPRKStepSetInterpolantDegree()`.  The user may request *k* in the
   range {0,..., *min(degree, kmax)*} where *kmax* depends on the choice of
   interpolation module. For Hermite interpolants *kmax = 5* and for Lagrange
   interpolants *kmax = 3*.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param t: the value of the independent variable at which the
        derivative is to be evaluated.
   :param k: the derivative order requested.
   :param dky: output vector (must be allocated by the user).

   :retval ARK_SUCCESS: if successful
   :retval ARK_BAD_K: if *k* is not in the range {0,..., *min(degree, kmax)*}.
   :retval ARK_BAD_T: if *t* is not in the interval :math:`[t_n-h_n, t_n]`
   :retval ARK_BAD_DKY: if the *dky* vector was ``NULL``
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``

   .. note::

      Dense outputs may or may not conserve the Hamiltonian. Our testing has
      shown that Lagrange interpolation typically performs well in this regard,
      while Hermite interpolation does not.

   .. warning::

      It is only legal to call this function after a successful
      return from :c:func:`SPRKStepEvolve()`.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetDky` instead.


.. _ARKODE.Usage.SPRKStep.OptionalOutputs:

Optional output functions
------------------------------


.. _ARKODE.Usage.SPRKStep.SPRKStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int SPRKStepGetNumSteps(void* arkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nsteps: number of steps taken in the solver.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumSteps` instead.


.. c:function:: int SPRKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param hlast: step size taken on the last internal step.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetLastStep` instead.


.. c:function:: int SPRKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)

   Returns the integration step size to be attempted on the next internal step.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param hcur: step size to be attempted on the next internal step.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetCurrentStep` instead.


.. c:function:: int SPRKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)

   Returns the current internal time reached by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tcur: current internal time reached.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetCurrentTime` instead.


.. c:function:: int SPRKStepGetCurrentState(void *arkode_mem, N_Vector *ycur)

   Returns the current internal solution reached by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ycur: current internal solution

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. warning::

      Users should exercise extreme caution when using this function,
      as altering values of *ycur* may lead to undesirable behavior, depending
      on the particular use case and on when this routine is called.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetCurrentState` instead.


.. c:function:: int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused, sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)

   Returns many of the most useful optional outputs in a single call.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nsteps: number of steps taken in the solver.
   :param hinused: actual value of initial step size.
   :param hlast: step size taken on the last internal step.
   :param hcur: step size to be attempted on the next internal step.
   :param tcur: current internal time reached.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetStepStats` instead.


.. c:function:: int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)

   Outputs all of the integrator statistics.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param outfile: pointer to output file.
   :param fmt: the output format:

       * :c:enumerator:`SUN_OUTPUTFORMAT_TABLE` -- prints a table of values
       * :c:enumerator:`SUN_OUTPUTFORMAT_CSV` -- prints a comma-separated list
         of key and value pairs e.g., ``key1,value1,key2,value2,...``

   :retval ARK_SUCCESS: -- if the output was successfully.
   :retval ARK_MEM_NULL: -- if the SPRKStep memory was ``NULL``.
   :retval ARK_ILL_INPUT: -- if an invalid formatting option was provided.

   .. note::

      The Python module ``tools/suntools`` provides utilities to read and output
      the data from a SUNDIALS CSV output file using the key and value pair
      format.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodePrintAllStats` instead.


.. c:function:: char *SPRKStepGetReturnFlagName(long int flag)

   Returns the name of the SPRKStep constant corresponding to *flag*.
   See :ref:`ARKODE.Constants`.

   :param flag: a return flag from an SPRKStep function.

   :returns: The return value is a string containing the name of the
             corresponding constant.

   .. warning::

      The user is responsible for freeing the returned string.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetReturnFlagName` instead.


.. c:function:: int SPRKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param step_attempts: number of steps attempted by solver.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumStepAttempts` instead.


.. c:function:: int SPRKStepGetNumRhsEvals(void* arkode_mem, long int* nf1, long int* nf2)

   Returns the number of calls to the user's right-hand
   side functions, :math:`f_1` and :math:`f_2` (so far).

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param nf1: number of calls to the user's :math:`f_1(t,p)` function.
   :param nf2: number of calls to the user's :math:`f_2(t,q)` function.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.2.0

      Use :c:func:`ARKodeGetNumRhsEvals` instead.


.. c:function:: int SPRKStepGetCurrentMethod(void* arkode_mem, ARKodeSPRKTable *sprk_table)

   Returns the SPRK method coefficient table currently in use by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param sprk_table: pointer to the SPRK method table.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``


.. c:function:: int SPRKStepGetUserData(void* arkode_mem, void** user_data)

   Returns the user data pointer previously set with
   :c:func:`SPRKStepSetUserData`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param user_data: memory reference to a user data pointer

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetUserData` instead.


.. _ARKODE.Usage.SPRKStep.SPRKStepRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int SPRKStepGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to have a root.

   For the components of :math:`g_i` for which a root was found, the sign of
   ``rootsfound[i]`` indicates the direction of zero-crossing. A value of +1
   indicates that :math:`g_i` is increasing, while a value of -1 indicates a
   decreasing :math:`g_i`.

   The user must allocate space for *rootsfound* prior to calling this function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param rootsfound: array of length *nrtfn* with the indices of the
        user functions :math:`g_i` found to have a root (the value of
        *nrtfn* was supplied in the call to
        :c:func:`SPRKStepRootInit()`).  For :math:`i = 0 \ldots`
        *nrtfn*-1, ``rootsfound[i]`` is nonzero if :math:`g_i` has a
        root, and 0 if not.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetRootInfo` instead.


.. c:function:: int SPRKStepGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param ngevals: number of calls made to :math:`g` so far.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumGEvals` instead.


.. _ARKODE.Usage.SPRKStep.SPRKStepExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int SPRKStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all SPRKStep solver parameters to the provided file pointer.

   The *fp* argument can be ``stdout`` or ``stderr``, or it may point to a
   specific file created using ``fopen``.

   When run in parallel, only one process should set a non-NULL value for this
   pointer, since parameters for all processes would be identical.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param fp: pointer to use for printing the solver parameters.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeWriteParameters` instead.


.. _ARKODE.Usage.SPRKStep.Reinitialization:

SPRKStep re-initialization function
-------------------------------------

To reinitialize the SPRKStep module for the solution of a new problem,
where a prior call to :c:func:`SPRKStepCreate` has been made, the
user must call the function :c:func:`SPRKStepReInit()`.  The new
problem must have the same size as the previous one.  This routine
retains the current settings for all SPRKStep module options and
performs the same input checking and initializations that are done in
:c:func:`SPRKStepCreate`, but it performs no memory allocation as it
assumes that the existing internal memory is sufficient for the new
problem.  A call to this re-initialization routine deletes the
solution history that was stored internally during the previous
integration, and deletes any previously-set *tstop* value specified via a
call to :c:func:`SPRKStepSetStopTime()`.  Following a successful call to
:c:func:`SPRKStepReInit()`, call :c:func:`SPRKStepEvolve()` again for the
solution of the new problem.

The use of :c:func:`SPRKStepReInit()` requires that the number of
Runge--Kutta stages, denoted by *s*, be no larger for the new problem than
for the previous problem.  This condition is automatically fulfilled
if the method order *q* is left unchanged.

One potential use of the :c:func:`SPRKStepReInit()` function is in the
treating of jump discontinuities in the RHS function :cite:p:`Tao:22`.
In lieu of including if statements within the RHS function to handle
discontinuities, it may be more computationally efficient to stop at each
point of discontinuity (e.g., through use of tstop or the rootfinding feature)
and restart the integrator with a readjusted ODE model, using a call to
this routine. We note that for the solution to retain temporal accuracy,
the RHS function should not incorporate the discontinuity.


.. c:function:: int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the SPRKStep
   time-stepper module.

   All previously set options are retained but may be updated by calling the
   appropriate "Set" functions.

   If an error occurred, :c:func:`SPRKStepReInit()` also sends an error message
   to the error handler function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param f1: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f1(t,q) = \frac{\partial V(t,q)}{\partial q}`
   :param f2: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f2(t,p) = \frac{\partial T(t,p)}{\partial p}`
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUT: if an argument had an illegal value.


.. _ARKODE.Usage.SPRKStep.Reset:

SPRKStep reset function
-----------------------

.. c:function:: int SPRKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)

   Resets the current SPRKStep time-stepper module state to the provided
   independent variable value and dependent variable vector.

   All previously set options are retained but may be updated by calling
   the appropriate "Set" functions.

   If an error occurred, :c:func:`SPRKStepReset()` also sends an error message
   to the error handler function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param tR: the value of the independent variable :math:`t`.
   :param yR: the value of the dependent variable vector :math:`y(t_R)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUTL: if an argument had an illegal value.

   .. note::

      By default the next call to :c:func:`SPRKStepEvolve()` will use the step
      size computed by SPRKStep prior to calling :c:func:`SPRKStepReset()`.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeReset` instead.
