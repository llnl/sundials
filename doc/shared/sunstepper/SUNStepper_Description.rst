.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @LLNL
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

.. _SUNStepper.Description:

The SUNStepper API
==================

.. versionadded:: 7.2.0

As with other SUNDIALS classes, the :c:type:`SUNStepper` abstract base class is
implemented using a C structure containing a ``content`` pointer to the derived
class member data and a structure of function pointers to the derived class
implementations of the virtual methods.

.. c:type:: SUNStepper

   An object for solving the IVP :eq:`SUNStepper_IVP`.

   The actual definition of the ``SUNStepper`` structure is kept private to
   allow for the object internals to change without impacting user code. The
   following sections describe the base class methods and the virtual methods
   that a must be provided by a derived class.

.. _SUNStepper.Description.BaseMethods:

Base Class Methods
------------------

This section describes methods provided by the :c:type:`SUNStepper` abstract
base class that aid the user in implementing derived classes. This includes
functions for creating and destroying a generic base class object, attaching and
retrieving the derived class ``content`` pointer, and setting function pointers
to derived class method implementations.

.. _SUNStepper.Description.BaseMethods.CreateDestroy:

Creating and Destroying an Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to creating an empty :c:type:`SUNStepper` using
:c:func:`SUNStepper_Create` described below, there is the
:c:func:`ARKodeCreateSUNStepper` function to construct a :c:type:`SUNStepper`
from an ARKODE integrator.

.. c:function:: SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper *stepper)

   This function creates a :c:type:`SUNStepper` object to which a user should
   attach the member data (content) pointer and method function pointers.

   :param sunctx: the SUNDIALS simulation context.
   :param stepper: a pointer to a stepper object.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   **Example usage:**

   .. code-block:: C

      /* create an instance of the base class */
      SUNStepper stepper = NULL;
      SUNErrCode err = SUNStepper_Create(sunctx, &stepper);

   .. note::

      See :numref:`SUNStepper.Description.BaseMethods.Content` and
      :numref:`SUNStepper.Description.BaseMethods.AttachFunctions`
      for details on how to attach member data and method function pointers.


.. c:function:: SUNErrCode SUNStepper_Destroy(SUNStepper *stepper)

   This function frees memory allocated by the :c:type:`SUNStepper` base class
   and uses the function pointer optionally specified with
   :c:func:`SUNStepper_SetDestroyFn` to free the content.

   :param stepper: a pointer to a stepper object.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      This function only frees memory allocated within the base class and the
      base class structure itself. The user is responsible for freeing any
      memory allocated for the member data (content).


.. _SUNStepper.Description.BaseMethods.SteppingFunctions:

Stepping Functions
^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_Evolve(SUNStepper stepper, sunrealtype tout, N_Vector vret, sunrealtype* tret)

   This function evolves the ODE :eq:`SUNStepper_IVP` towards the time ``tout``
   and stores the solution at time ``tret`` in ``vret``.

   :param stepper: the stepper object.
   :param tout: the time to evolve towards.
   :param vret: on output, the state at time ``tret``.
   :param tret: the time corresponding to the output value ``vret``.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNStepper_OneStep(SUNStepper stepper, sunrealtype tout, N_Vector vret, sunrealtype* tret)

   This function evolves the ODE :eq:`SUNStepper_IVP` *one timestep* towards
   the time ``tout`` and stores the solution at time ``tret`` in ``vret``.

   :param stepper: the stepper object.
   :param tout: the time to evolve towards.
   :param vret: on output, the state at time ``tret``.
   :param tret: the time corresponding to the output value ``vret``.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNStepper_FullRhs(SUNStepper stepper, sunrealtype t, N_Vector v, N_Vector f, SUNFullRhsMode mode)

   This function computes the full right-hand side function of the ODE,
   :math:`f(t, v) + r(t)` in :eq:`SUNStepper_IVP` for a given value of the
   independent variable ``t`` and state vector ``v``.

   :param stepper: the stepper object.
   :param t: the current value of the independent variable.
   :param v: the current value of the dependent variable vector.
   :param f: the output vector for the ODE right-hand side,
      :math:`f(t, v) + r(t)`, in :eq:`SUNStepper_IVP`.
   :param mode: the purpose of the right-hand side evaluation.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_ReInit(SUNStepper stepper, sunrealtype t0, N_Vector v0)

   This function reinitalizes the stepper to solve a new problem with the given initial
   condition and clears all counters.

   :param stepper: the stepper object.
   :param t0: the value of the independent variable :math:`t_0`.
   :param v0: the value of the dependent variable vector :math:`v(t_0)`.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 7.3.0


.. c:function:: SUNErrCode SUNStepper_Reset(SUNStepper stepper, sunrealtype tR, N_Vector vR)

   This function resets the stepper state to the provided independent variable
   value and dependent variable vector.

   :param stepper: the stepper object.
   :param tR: the value of the independent variable :math:`t_R`.
   :param vR: the value of the dependent variable vector :math:`v(t_R)`.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_ResetCheckpointIndex(SUNStepper stepper, suncountertype ckptIdxR)

   This function resets the index at which new checkpoints will be inserted to `ckptIdxR`.

   :param stepper: the stepper object.
   :param ckptIdxR: the step index to begin checkpointing from
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 7.3.0


.. c:function:: SUNErrCode SUNStepper_SetStopTime(SUNStepper stepper, sunrealtype tstop)

   This function specifies the value of the independent variable :math:`t` past
   which the solution is not to proceed.

   :param stepper: the stepper object.
   :param tstop: stopping time for the stepper.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetStepDirection(SUNStepper stepper, sunrealtype stepdir)

   This function specifies the direction of integration (forward or backward).

   :param stepper: the stepper object.
   :param stepdir: value whose sign determines the direction. A positive value
      selects forward integration, a negative value selects backward
      integration, and zero leaves the current direction unchanged.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetForcing(SUNStepper stepper, sunrealtype tshift, sunrealtype tscale, N_Vector* forcing, int nforcing)

   This function sets the data necessary to compute the forcing term
   :eq:`SUNStepper_forcing`. This includes the shift and scaling factors for the
   normalized time :math:`\frac{t - t_{\text{shift}}}{t_{\text{scale}}}` and the
   array of polynomial coefficient vectors :math:`\widehat{f}_k`.

   :param stepper: a stepper object.
   :param tshift: the time shift to apply to the current time when computing
      the forcing, :math:`t_{\text{shift}}`.
   :param tscale: the time scaling to apply to the current time when computing
      the forcing, :math:`t_{\text{scale}}`.
   :param forcing: a pointer to an array of forcing vectors,
      :math:`\widehat{f}_k`.
   :param nforcing: the number of forcing vectors, :math:`n_{\text{forcing}}`. A
      value of 0 effectively eliminates the forcing term.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      When integrating the ODE :eq:`SUNStepper_IVP` the :c:type:`SUNStepper` is
      responsible for evaluating ODE right-hand side function :math:`f(t, v)` as
      well as computing and applying the forcing term :eq:`SUNStepper_forcing`
      to obtain the full right-hand side of the ODE :eq:`SUNStepper_IVP`.

.. c:function:: SUNErrCode SUNStepper_GetNumSteps(SUNStepper stepper, suncountertype* nst)

   This function gets the number of successful time steps taken by the stepper
   since it was last initialized.

   :param stepper: the stepper object.
   :param nst: on output, the number of time steps.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 7.3.0


.. _SUNStepper.Description.BaseMethods.RhsMode:

The Right-Hand Side Evaluation Mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:enum:: SUNFullRhsMode

   A flag indicating the purpose of a right-hand side function evaluation.

   .. c:enumerator:: SUN_FULLRHS_START

      Evaluate at the beginning of the simulation.

   .. c:enumerator:: SUN_FULLRHS_END

      Evaluate at the end of a successful step.

   .. c:enumerator:: SUN_FULLRHS_OTHER

      Evaluate elsewhere, e.g., for dense output.


.. _SUNStepper.Description.BaseMethods.Content:

Attaching and Accessing the Content Pointer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void *content)

   This function attaches a member data (content) pointer to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param content: a pointer to the stepper member data.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void **content)

   This function retrieves the member data (content) pointer from a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param content: a pointer to set to the stepper member data pointer.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


Handling Warnings and Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An implementation of a :c:type:`SUNStepper` may have a system of warning and error
handling that cannot be encoded as a :c:type:`SUNErrCode` which is the return
type of all :c:type:`SUNStepper` functions. Therefore, we provide the following
function to get and set a separate flag associated with a stepper.

.. c:function:: SUNErrCode SUNStepper_SetLastFlag(SUNStepper stepper, int last_flag)

   This function sets a flag that can be used by :c:type:`SUNStepper` implementations to
   indicate warnings or errors that occurred during an operation, e.g.,
   :c:func:`SUNStepper_Evolve`.

   :param stepper: the stepper object.
   :param last_flag: the flag value.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNStepper_GetLastFlag(SUNStepper stepper, int *last_flag)

   This function provides the last value of the flag used by the :c:type:`SUNStepper`
   implementation to indicate warnings or errors that occurred during an
   operation, e.g., :c:func:`SUNStepper_Evolve`.

   :param stepper: the stepper object.
   :param last_flag: A pointer to where the flag value will be written.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. _SUNStepper.Description.BaseMethods.AttachFunctions:

Setting Member Functions
^^^^^^^^^^^^^^^^^^^^^^^^

The functions in this section are used to specify how each operation on a
:c:type:`SUNStepper` implementation is performed. Technically, all of these
functions are optional to call; the functions that need to be attached are
determined by the "consumer" of the :c:type:`SUNStepper`.

.. c:function:: SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn)

   This function attaches a :c:type:`SUNStepperEvolveFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperEvolveFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetOneStepFn(SUNStepper stepper, SUNStepperOneStepFn fn)

   This function attaches a :c:type:`SUNStepperOneStepFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperOneStepFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn)

   This function attaches a :c:type:`SUNStepperFullRhsFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperFullRhsFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNStepper_SetReInitFn(SUNStepper stepper, SUNStepperResetFn fn)

   This function attaches a :c:type:`SUNStepperReInitFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperReInitFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 7.3.0

.. c:function:: SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn)

   This function attaches a :c:type:`SUNStepperResetFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperResetFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNStepper_SetResetCheckpointIndexFn(SUNStepper stepper, SUNStepperResetCheckpointIndexFn fn)

   This function attaches a :c:type:`SUNStepperResetCheckpointIndexFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperResetCheckpointIndexFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 7.3.0


.. c:function:: SUNErrCode SUNStepper_SetStopTimeFn(SUNStepper stepper, SUNStepperSetStopTimeFn fn)

   This function attaches a :c:type:`SUNStepperSetStopTimeFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperSetStopTimeFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetStepDirectionFn(SUNStepper stepper, SUNStepperSetStepDirectionFn fn)

   This function attaches a :c:type:`SUNStepperSetStepDirectionFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperSetStepDirectionFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetForcingFn(SUNStepper stepper, SUNStepperSetForcingFn fn)

   This function attaches a :c:type:`SUNStepperSetForcingFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperSetForcingFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNStepper_SetGetNumStepsFn(SUNStepper stepper, SUNStepperGetNumStepsFn fn)

   This function attaches a :c:type:`SUNStepperGetNumStepsFn` function to a
   :c:type:`SUNStepper` object.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperGetNumStepsFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 7.3.0


.. c:function:: SUNErrCode SUNStepper_SetDestroyFn(SUNStepper stepper, SUNStepperDestroyFn fn)

   This function attaches a :c:type:`SUNStepperDestroyFn` function to a
   :c:type:`SUNStepper`. The provided function is responsible for freeing any
   memory allocated for the :c:type:`SUNStepper` content.

   :param stepper: a stepper object.
   :param fn: the :c:type:`SUNStepperDestroyFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. _SUNStepper.Description.ImplMethods:

Implementation Specific Methods
-------------------------------

This section describes the virtual methods defined by the :c:type:`SUNStepper`
abstract base class.


.. c:type:: SUNErrCode (*SUNStepperEvolveFn)(SUNStepper stepper, sunrealtype tout, N_Vector vret, sunrealtype* tret)

   This type represents a function with the signature of
   :c:func:`SUNStepper_Evolve`.


.. c:type:: SUNErrCode (*SUNStepperOneStepFn)(SUNStepper stepper, sunrealtype tout, N_Vector vret, sunrealtype* tret)

   This type represents a function with the signature of
   :c:func:`SUNStepper_OneStep`.


.. c:type:: SUNErrCode (*SUNStepperFullRhsFn)(SUNStepper stepper, sunrealtype t, N_Vector v, N_Vector f, SUNFullRhsMode mode)

   This type represents a function with the signature of
   :c:func:`SUNStepper_FullRhs`.


.. c:type:: SUNErrCode (*SUNStepperResetFn)(SUNStepper stepper, sunrealtype tR, N_Vector vR)

   This type represents a function with the signature of
   :c:func:`SUNStepper_Reset`.


.. c:type:: SUNErrCode (*SUNStepperReInitFn)(SUNStepper stepper, sunrealtype tR, N_Vector vR)

   This type represents a function with the signature of
   :c:func:`SUNStepper_ReInit`.

   .. versionadded:: 7.3.0


.. c:type:: SUNErrCode (*SUNStepperResetCheckpointIndexFn)(SUNStepper stepper, suncountertype ckptIdxR)

   This type represents a function with the signature of
   :c:func:`SUNStepper_ResetCheckpointIndex`.

   .. versionadded:: 7.3.0


.. c:type:: SUNErrCode (*SUNStepperSetStopTimeFn)(SUNStepper stepper, sunrealtype tstop)

   This type represents a function with the signature of
   :c:func:`SUNStepper_SetStopTime`.


.. c:type:: SUNErrCode (*SUNStepperSetStepDirectionFn)(SUNStepper stepper, sunrealtype stepdir)

   This type represents a function with the signature of
   :c:func:`SUNStepper_SetStepDirection`.


.. c:type:: SUNErrCode (*SUNStepperSetForcingFn)(SUNStepper stepper, sunrealtype tshift, sunrealtype tscale, N_Vector* forcing, int nforcing)

   This type represents a function with the signature of
   :c:func:`SUNStepper_SetForcing`.

.. c:type:: SUNErrCode (*SUNStepperDestroyFn)(SUNStepper stepper)

   This type represents a function with the signature similar to
   :c:func:`SUNStepper_Destroy` for freeing the content associated with a
   :c:type:`SUNStepper`.

.. c:type:: SUNErrCode (*SUNStepperGetNumStepsFn)(SUNStepper stepper, suncountertype* nst)

   This type represents a function with the signature of
   :c:func:`SUNStepper_GetNumSteps`.

   .. versionadded:: 7.3.0
