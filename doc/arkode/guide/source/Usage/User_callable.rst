.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ UMBC
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

.. _ARKODE.Usage.UserCallable:

ARKODE User-callable functions
================================

This section describes the shared ARKODE functions that are called by
the user to setup and then solve an IVP. Some of these are required;
however, starting with :numref:`ARKODE.Usage.OptionalInputs`,
the functions listed involve optional inputs/outputs or restarting,
and those paragraphs may be skipped for a casual use of ARKODE.
In any case, refer to the preceding section,
:numref:`ARKODE.Usage.Skeleton`, for the correct order of these calls.

On an error, each user-callable function returns a negative value (or
``NULL`` if the function returns a pointer) and sends an error message
to the error handler, which prints the message to ``stderr`` by default.
However, the user can set a file as error output or can
provide their own error handler (see :numref:`SUNDIALS.Errors` for details).

We note that depending on the choice of time-stepping module, only a
subset of ARKODE's user-callable functions will be applicable/supported.
We thus categorize the functions below into five groups:

A. functions that apply for all time-stepping modules,

B. functions that apply for time-stepping modules that allow temporal adaptivity,

C. functions that apply for time-stepping modules that utilize implicit solvers (nonlinear or linear),

D. functions that apply for time-stepping modules that support non-identity mass matrices, and

E. functions that apply for time-stepping modules that support relaxation Runge--Kutta methods.

In the function descriptions below, we identify those that have any of the restrictions B-E above.
Then in the introduction for each of the stepper-specific documentation sections
(:numref:`ARKODE.Usage.ARKStep.UserCallable`,
:numref:`ARKODE.Usage.ERKStep.UserCallable`,
:numref:`ARKODE.Usage.ForcingStep.UserCallable`,
:numref:`ARKODE.Usage.LSRKStep.UserCallable`,
:numref:`ARKODE.Usage.MRIStep.UserCallable`,
:numref:`ARKODE.Usage.SplittingStep.UserCallable`,
and :numref:`ARKODE.Usage.SPRKStep.UserCallable`)
we clarify the categories of these functions that are supported.


.. _ARKODE.Usage.Initialization:

ARKODE initialization and deallocation functions
------------------------------------------------------

For functions to create an ARKODE stepper instance see :c:func:`ARKStepCreate`,
:c:func:`ERKStepCreate`, :c:func:`ForcingStepCreate`,
:c:func:`LSRKStepCreateSTS`, :c:func:`LSRKStepCreateSSP`,
:c:func:`MRIStepCreate`, :c:func:`SplittingStepCreate`, or
:c:func:`SPRKStepCreate`.

.. c:function:: void ARKodeFree(void** arkode_mem)

   This function frees the problem memory created a stepper constructor.

   :param arkode_mem: pointer to the ARKODE stepper memory block.
   :return: none

   .. versionadded:: 6.1.0

      This function replaces stepper specific versions in ARKStep, ERKStep,
      MRIStep, and SPRKStep.


.. _ARKODE.Usage.Tolerances:

ARKODE tolerance specification functions
------------------------------------------------------

These functions specify the integration tolerances. One of them
**should** be called before the first call to
:c:func:`ARKodeEvolve`; otherwise default values of ``reltol =
1e-4`` and ``abstol = 1e-9`` will be used, which may be entirely
incorrect for a specific problem.

The integration tolerances ``reltol`` and ``abstol`` define a vector
of error weights, ``ewt``.  In the case of
:c:func:`ARKodeSStolerances`, this vector has components

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol);

whereas in the case of :c:func:`ARKodeSVtolerances` the vector components
are given by

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol[i]);

This vector is used in all error and convergence tests, which use a
weighted RMS norm on all error-like vectors :math:`v`:

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N (v_i\; ewt_i)^2 \right)^{1/2},

where :math:`N` is the problem dimension.

Alternatively, the user may supply a custom function to supply the
``ewt`` vector, through a call to :c:func:`ARKodeWFtolerances`.



.. c:function:: int ARKodeSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)

   This function specifies scalar relative and absolute tolerances.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param reltol: scalar relative tolerance.
   :param abstol: scalar absolute tolerance.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC:  ``arkode_mem`` was not allocated.
   :retval ARK_ILL_INPUT: an argument had an illegal value (e.g. a negative tolerance).

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.scalar_tolerances".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)

   This function specifies a scalar relative tolerance and a vector
   absolute tolerance (a potentially different absolute tolerance for
   each vector component).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param reltol: scalar relative tolerance.
   :param abstol: vector containing the absolute tolerances for each
                  solution component.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC:  ``arkode_mem`` was not allocated.
   :retval ARK_ILL_INPUT: an argument had an illegal value (e.g. a negative tolerance).

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeWFtolerances(void* arkode_mem, ARKEwtFn efun)

   This function specifies a user-supplied function *efun* to compute
   the error weight vector ``ewt``.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param efun: the name of the function (of type :c:func:`ARKEwtFn`)
                that implements the error weight vector computation.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC:  ``arkode_mem`` was not allocated.

   .. versionadded:: 6.1.0


Moreover, for problems involving a non-identity mass matrix
:math:`M \ne I`, the units of the solution vector :math:`y` may differ
from the units of the IVP, posed for the vector :math:`My`.  When this
occurs, iterative solvers for the Newton linear systems and the mass
matrix linear systems may require a different set of tolerances.
Since the relative tolerance is dimensionless, but the absolute
tolerance encodes a measure of what is "small" in the units of the
respective quantity, a user may optionally define absolute tolerances
in the equation units.  In this case, ARKODE defines a vector of residual
weights, ``rwt`` for measuring convergence of these iterative solvers.
In the case of :c:func:`ARKodeResStolerance`, this vector has components

.. code-block:: c

   rwt[i] = 1.0/(reltol*abs(My[i]) + rabstol);

whereas in the case of :c:func:`ARKodeResVtolerance` the vector components
are given by

.. code-block:: c

   rwt[i] = 1.0/(reltol*abs(My[i]) + rabstol[i]);

This residual weight vector is used in all iterative solver
convergence tests, which similarly use a weighted RMS norm on all
residual-like vectors :math:`v`:

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N (v_i\; rwt_i)^2 \right)^{1/2},

where :math:`N` is the problem dimension.

As with the error weight vector, the user may supply a custom function
to supply the ``rwt`` vector, through a call to
:c:func:`ARKodeResFtolerance`.  Further information on all three of
these functions is provided below.



.. c:function:: int ARKodeResStolerance(void* arkode_mem, sunrealtype rabstol)

   This function specifies a scalar absolute residual tolerance.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rabstol: scalar absolute residual tolerance.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC:  ``arkode_mem`` was not allocated.
   :retval ARK_ILL_INPUT: an argument had an illegal value (e.g. a negative tolerance).

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeResVtolerance(void* arkode_mem, N_Vector rabstol)

   This function specifies a vector of absolute residual tolerances.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rabstol: vector containing the absolute residual
                   tolerances for each solution component.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC:  ``arkode_mem`` was not allocated.
   :retval ARK_ILL_INPUT: an argument had an illegal value (e.g. a negative tolerance).

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeResFtolerance(void* arkode_mem, ARKRwtFn rfun)

   This function specifies a user-supplied function *rfun* to compute
   the residual weight vector ``rwt``.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rfun: the name of the function (of type :c:func:`ARKRwtFn`)
                that implements the residual weight vector computation.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC:  ``arkode_mem`` was not allocated.

   .. versionadded:: 6.1.0


General advice on the choice of tolerances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many users, the appropriate choices for tolerance values in
``reltol``, ``abstol``, and ``rabstol`` are a concern. The following pieces
of advice are relevant.

(1) The scalar relative tolerance ``reltol`` is to be set to control
    relative errors. So a value of :math:`10^{-4}` means that errors
    are controlled to .01%. We do not recommend using ``reltol`` larger
    than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so
    small that it is comparable to the unit roundoff of the machine
    arithmetic (generally around :math:`10^{-15}` for double-precision).

(2) The absolute tolerances ``abstol`` (whether scalar or vector) need
    to be set to control absolute errors when any components of the
    solution vector :math:`y` may be so small that pure relative error
    control is meaningless.  For example, if :math:`y_i` starts at some
    nonzero value, but in time decays to zero, then pure relative
    error control on :math:`y_i` makes no sense (and is overly costly)
    after :math:`y_i` is below some noise level. Then ``abstol`` (if
    scalar) or ``abstol[i]`` (if a vector) needs to be set to that
    noise level. If the different components have different noise
    levels, then ``abstol`` should be a vector.  For example, see the
    example problem ``ark_robertson.c``, and the discussion
    of it in the ARKODE Examples Documentation :cite:p:`arkode_ex`.  In that
    problem, the three components vary between 0 and 1, and have
    different noise levels; hence the ``atols`` vector therein. It is
    impossible to give any general advice on ``abstol`` values,
    because the appropriate noise levels are completely
    problem-dependent. The user or modeler hopefully has some idea as
    to what those noise levels are.

(3) The residual absolute tolerances ``rabstol`` (whether scalar or
    vector) follow a similar explanation as for ``abstol``, except
    that these should be set to the noise level of the equation
    components, i.e. the noise level of :math:`My`.  For problems in
    which :math:`M=I`, it is recommended that ``rabstol`` be left
    unset, which will default to the already-supplied ``abstol``
    values.

(4) Finally, it is important to pick all the tolerance values
    conservatively, because they control the error committed on each
    individual step. The final (global) errors are an accumulation of
    those per-step errors, where that accumulation factor is
    problem-dependent.  A general rule of thumb is to reduce the
    tolerances by a factor of 10 from the actual desired limits on
    errors.  So if you want .01% relative accuracy (globally), a good
    choice for ``reltol`` is :math:`10^{-5}`.  In any case, it is
    a good idea to do a few experiments with the tolerances to see how
    the computed solution values vary as tolerances are reduced.



Advice on controlling nonphysical negative values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications, some components in the true solution are always
positive or non-negative, though at times very small.  In the
numerical solution, however, small negative (nonphysical) values
can then occur. In most cases, these values are harmless, and simply
need to be controlled, not eliminated, but in other cases any value
that violates a constraint may cause a simulation to halt. For both of
these scenarios the following pieces of advice are relevant.

(1) The best way to control the size of unwanted negative computed
    values is with tighter absolute tolerances.  Again this requires
    some knowledge of the noise level of these components, which may
    or may not be different for different components. Some
    experimentation may be needed.

(2) If output plots or tables are being generated, and it is important
    to avoid having negative numbers appear there (for the sake of
    avoiding a long explanation of them, if nothing else), then
    eliminate them, but only in the context of the output medium. Then
    the internal values carried by the solver are unaffected. Remember
    that a small negative value in :math:`y` returned by ARKODE, with
    magnitude comparable to ``abstol`` or less, is equivalent to zero
    as far as the computation is concerned.

(3) The user's right-hand side routines :math:`f^E` and :math:`f^I`
    should never change a negative value in the solution vector :math:`y`
    to a non-negative value in attempt to "fix" this problem,
    since this can lead to numerical instability.  If the :math:`f^E`
    or :math:`f^I` routines cannot tolerate a zero or negative value
    (e.g. because there is a square root or log), then the offending
    value should be changed to zero or a tiny positive number in a
    temporary variable (not in the input :math:`y` vector) for the
    purposes of computing :math:`f^E(t, y)` or :math:`f^I(t, y)`.

(4) Some of ARKODE's time stepping modules support component-wise
    constraints on solution components, :math:`y_i < 0`,
    :math:`y_i \le 0`, :math:`y_i > 0`, or :math:`y_i \ge 0`, through
    the user-callable function :c:func:`ARKodeSetConstraints`.  At each
    internal time step, if any constraint is violated then ARKODE will
    attempt a smaller time step that should not violate this constraint.
    This reduced step size is chosen such that the step size is the
    largest possible but where the solution component satisfies the
    constraint.

(5) For time-stepping modules that support temporal adaptivity,
    positivity and non-negativity constraints on components can also be
    enforced by use of the recoverable error return feature in the
    user-supplied right-hand side function(s). When a recoverable error
    is encountered, ARKODE will retry the step with a smaller step size,
    which typically alleviates the problem.  However, since this reduced
    step size is chosen without knowledge of the solution constraint, it
    may be overly conservative.  Thus this option involves some additional
    overhead cost, and should only be exercised if the above recommendations
    are unsuccessful.



.. _ARKODE.Usage.LinearSolvers:

Linear solver interface functions
-------------------------------------------

As previously explained, the Newton iterations used in solving
implicit systems within ARKODE require the solution of linear
systems of the form

.. math::
   \mathcal{A}\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right)

where

.. math::
   \mathcal{A} \approx M - \gamma J, \qquad J = \frac{\partial f^I}{\partial y}.

ARKODE's ARKLS linear solver interface supports all valid
``SUNLinearSolver`` modules for this task.

Matrix-based ``SUNLinearSolver`` modules utilize ``SUNMatrix`` objects
to store the approximate Jacobian matrix :math:`J`, the Newton matrix
:math:`\mathcal{A}`, the mass matrix :math:`M`, and, when using direct
solvers, the factorizations used throughout the solution process.

Matrix-free ``SUNLinearSolver`` modules instead use iterative methods
to solve the Newton systems of equations, and only require the
*action* of the matrix on a vector, :math:`\mathcal{A}v`.  With most
of these methods, preconditioning can be done on the left only, on the
right only, on both the left and the right, or not at all.  The
exceptions to this rule are SPFGMR that supports right preconditioning
only and PCG that performs symmetric preconditioning.  For the
specification of a preconditioner, see the iterative linear solver
portions of :numref:`ARKODE.Usage.OptionalInputs` and
:numref:`ARKODE.Usage.UserSupplied`.

If preconditioning is done, user-supplied functions should be used to
define left and right preconditioner matrices :math:`P_1` and
:math:`P_2` (either of which could be the identity matrix), such that
the product :math:`P_{1}P_{2}` approximates the Newton matrix
:math:`\mathcal{A} = M - \gamma J`.

To specify a generic linear solver for ARKODE to use for the Newton
systems, after the call to ``*StepCreate`` but before any
calls to :c:func:`ARKodeEvolve`, the user's program must create the
appropriate ``SUNLinearSolver`` object and call the function
:c:func:`ARKodeSetLinearSolver`, as documented below.  To create
the ``SUNLinearSolver`` object, the user may call one of the
SUNDIALS-packaged SUNLinSol module constructor routines via a call of
the form

.. code:: c

   SUNLinearSolver LS = SUNLinSol_*(...);

The current list of SUNDIALS-packaged SUNLinSol modules, and their
constructor routines, may be found in chapter :numref:`SUNLinSol`.
Alternately, a user-supplied ``SUNLinearSolver`` module may be created
and used.  Specific information on how to create such user-provided
modules may be found in :numref:`SUNLinSol.API.Custom`.

Once this solver object has been constructed, the user should attach
it to ARKODE via a call to :c:func:`ARKodeSetLinearSolver`. The
first argument passed to this function is the ARKODE memory pointer
returned by ``*StepCreate``; the second argument is the
``SUNLinearSolver`` object created above.  The third argument is an
optional ``SUNMatrix`` object to accompany matrix-based
``SUNLinearSolver`` inputs (for matrix-free linear solvers, the third
argument should be ``NULL``).  A call to this function initializes the
ARKLS linear solver interface, linking it to the ARKODE integrator,
and allows the user to specify additional parameters and routines
pertinent to their choice of linear solver.

.. c:function:: int ARKodeSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix J)

   This function specifies the ``SUNLinearSolver`` object that ARKODE
   should use, as well as a template Jacobian ``SUNMatrix`` object (if
   applicable).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param LS: the ``SUNLinearSolver`` object to use.
   :param J: the template Jacobian ``SUNMatrix`` object to use (or
             ``NULL`` if not applicable).

   :retval ARKLS_SUCCESS:   the function exited successfully.
   :retval ARKLS_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARKLS_MEM_FAIL:  there was a memory allocation failure.
   :retval ARKLS_ILL_INPUT: ARKLS is incompatible with the
                            provided *LS* or *J* input objects, or the current
                            ``N_Vector`` module.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      If *LS* is a matrix-free linear solver, then the *J*
      argument should be ``NULL``.

      If *LS* is a matrix-based linear solver, then the template Jacobian
      matrix *J* will be used in the solve process, so if additional
      storage is required within the ``SUNMatrix`` object (e.g. for
      factorization of a banded matrix), ensure that the input object is
      allocated with sufficient size (see the documentation of
      the particular SUNMATRIX type in the :numref:`SUNMatrix` for
      further information).

      When using sparse linear solvers, it is typically much more
      efficient to supply *J* so that it includes the full sparsity
      pattern of the Newton system matrices :math:`\mathcal{A} =
      M-\gamma J`, even if *J* itself has zeros in nonzero
      locations of :math:`M`.  The reasoning for this is
      that :math:`\mathcal{A}` is constructed in-place, on top of the
      user-specified values of *J*, so if the sparsity pattern in *J* is
      insufficient to store :math:`\mathcal{A}` then it will need to be
      resized internally by ARKODE.

   .. versionadded:: 6.1.0





.. _ARKODE.Usage.MassMatrixSolvers:

Mass matrix solver specification functions
-------------------------------------------

As discussed in :numref:`ARKODE.Mathematics.MassSolve`, if the ODE
system involves a non-identity mass matrix :math:`M\ne I`, then ARKODE
must solve linear systems of the form

.. math::
    M x = b.

ARKODE's ARKLS mass-matrix linear solver interface supports all valid
``SUNLinearSolver`` modules for this task.  For iterative linear
solvers, user-supplied preconditioning can be applied.  For the
specification of a preconditioner, see the iterative linear solver
portions of :numref:`ARKODE.Usage.OptionalInputs` and
:numref:`ARKODE.Usage.UserSupplied`.  If preconditioning is to be
performed, user-supplied functions should be used to define left and
right preconditioner matrices :math:`P_1` and :math:`P_2` (either of
which could be the identity matrix), such that the product
:math:`P_{1}P_{2}` approximates the mass matrix :math:`M`.

To specify a generic linear solver for ARKODE to use for mass matrix
systems, after the call to ``*StepCreate`` but before any
calls to :c:func:`ARKodeEvolve`, the user's program must create the
appropriate ``SUNLinearSolver`` object and call the function
:c:func:`ARKodeSetMassLinearSolver`, as documented below.  The
first argument passed to this function is the ARKODE memory
pointer returned by ``*StepCreate``; the second argument is
the desired ``SUNLinearSolver`` object to use for solving mass matrix
systems.  The third object is a template ``SUNMatrix`` to use with the
provided ``SUNLinearSolver`` (if applicable).  The fourth input is a
flag to indicate whether the mass matrix is time-dependent,
i.e. :math:`M = M(t)`, or not.  A call to this function initializes the
ARKLS mass matrix linear solver interface, linking this to the main
ARKODE integrator, and allows the user to specify additional
parameters and routines pertinent to their choice of linear solver.

Note: if the user program includes linear solvers for *both* the
Newton and mass matrix systems, these must have the same type:

* If both are matrix-based, then they must utilize the same
  ``SUNMatrix`` type, since these will be added when forming the
  Newton system matrix :math:`\mathcal{A}`.  In this case, both the
  Newton and mass matrix linear solver interfaces can use the same
  ``SUNLinearSolver`` object, although different solver objects
  (e.g. with different solver parameters) are also allowed.

* If both are matrix-free, then the Newton and mass matrix
  ``SUNLinearSolver`` objects must be different.  These may even use
  different solver algorithms (SPGMR, SPBCGS, etc.), if desired.
  For example, if the mass matrix is symmetric but the Jacobian is not,
  then PCG may be used for the mass matrix systems and SPGMR for the
  Newton systems.


.. c:function:: int ARKodeSetMassLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix M, sunbooleantype time_dep)

   This function specifies the ``SUNLinearSolver`` object
   that ARKODE should use for mass matrix systems, as well as a
   template ``SUNMatrix`` object.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param LS: the ``SUNLinearSolver`` object to use.
   :param M: the template mass ``SUNMatrix`` object to use.
   :param time_dep: flag denoting whether the mass matrix depends on
                    the independent variable (:math:`M = M(t)`) or not (:math:`M
                    \ne M(t)`).  ``SUNTRUE`` indicates time-dependence of the
                    mass matrix.

   :retval ARKLS_SUCCESS:   the function exited successfully.
   :retval ARKLS_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARKLS_MEM_FAIL:  there was a memory allocation failure.
   :retval ARKLS_ILL_INPUT: ARKLS is incompatible with the
                            provided *LS* or *M* input objects, or the current
                            ``N_Vector`` module.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      If *LS* is a matrix-free linear solver, then the *M*
      argument should be ``NULL``.

      If *LS* is a matrix-based linear solver, then the template mass
      matrix *M* will be used in the solve process, so if additional
      storage is required within the ``SUNMatrix`` object (e.g. for
      factorization of a banded matrix), ensure that the input object is
      allocated with sufficient size.

      If called with *time_dep* set to ``SUNFALSE``, then the mass matrix is
      only computed and factored once (or when either ``*StepReInit``
      or :c:func:`ARKodeResize` are called), with the results reused
      throughout the entire ARKODE simulation.

      Unlike the system Jacobian, the system mass matrix is not approximated
      using finite-differences of any functions provided to ARKODE.  Hence,
      use of the a matrix-based *LS* requires the user to provide a
      mass-matrix constructor routine (see :c:type:`ARKLsMassFn` and
      :c:func:`ARKodeSetMassFn`).

      Similarly, the system mass matrix-vector-product is not approximated
      using finite-differences of any functions provided to ARKODE.  Hence,
      use of a matrix-free *LS* requires the user to provide a
      mass-matrix-times-vector product routine (see
      :c:type:`ARKLsMassTimesVecFn` and :c:func:`ARKodeSetMassTimes`).

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.NonlinearSolvers:

Nonlinear solver interface functions
-------------------------------------------

When changing the nonlinear solver in ARKODE, after the
call to ``*StepCreate`` but before any calls to
:c:func:`ARKodeEvolve`, the user's program must create the
appropriate ``SUNNonlinearSolver`` object and call
:c:func:`ARKodeSetNonlinearSolver`, as documented below.  If any
calls to :c:func:`ARKodeEvolve` have been made, then ARKODE will
need to be reinitialized by calling ``*StepReInit`` to
ensure that the nonlinear solver is initialized correctly before any
subsequent calls to :c:func:`ARKodeEvolve`.

The first argument passed to the routine
:c:func:`ARKodeSetNonlinearSolver` is the ARKODE memory pointer
returned by ``*StepCreate``; the second argument passed
to this function is the desired ``SUNNonlinearSolver`` object to use for
solving the nonlinear system for each implicit stage. A call to this
function attaches the nonlinear solver to the main ARKODE integrator.


.. c:function:: int ARKodeSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS)

   This function specifies the ``SUNNonlinearSolver`` object
   that ARKODE should use for implicit stage solves.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param NLS: the ``SUNNonlinearSolver`` object to use.

   :retval ARK_SUCCESS:   the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_MEM_FAIL:  there was a memory allocation failure.
   :retval ARK_ILL_INPUT: ARKODE is incompatible with the
                          provided *NLS* input object.
   :retval ARK_STEPPER_UNSUPPORTED: nonlinear solvers are not supported by
                                    the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      ARKODE will use the Newton ``SUNNonlinearSolver`` module by
      default; a call to this routine replaces that module with the
      supplied *NLS* object.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.RootFinding:

Rootfinding initialization function
--------------------------------------

As described in :numref:`ARKODE.Mathematics.Rootfinding`, while
solving the IVP, ARKODE's time-stepping modules have the capability to
find the roots of a set of user-defined functions.  To activate the
root-finding algorithm, call the following function.  This is normally
called only once, prior to the first call to
:c:func:`ARKodeEvolve`, but if the rootfinding problem is to be
changed during the solution, :c:func:`ARKodeRootInit` can also be
called prior to a continuation call to :c:func:`ARKodeEvolve`.

.. note::

   The solution is interpolated to the times at which roots are found.


.. c:function:: int ARKodeRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system.  It must be called after
   ``*StepCreate``, and before :c:func:`ARKodeEvolve`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nrtfn: number of functions :math:`g_i`, an integer :math:`\ge` 0.
   :param g: name of user-supplied function, of type :c:func:`ARKRootFn`,
             defining the functions :math:`g_i` whose roots are sought.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_MEM_FAIL:  there was a memory allocation failure.
   :retval ARK_ILL_INPUT: *nrtfn* is greater than zero but *g* is ``NULL``.

   .. note::

      To disable the rootfinding feature after it has already
      been initialized, or to free memory associated with ARKODE's
      rootfinding module, call *ARKodeRootInit* with *nrtfn = 0*.

      Similarly, if a new IVP is to be solved with a call to
      ``*StepReInit``, where the new IVP has no rootfinding
      problem but the prior one did, then call *ARKodeRootInit* with
      *nrtfn = 0*.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.Integration:

ARKODE solver function
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  The input argument *itask* specifies one of two
modes as to where ARKODE is to return a solution.  These modes are modified if
the user has set a stop time (with a call to the optional input function
:c:func:`ARKodeSetStopTime`) or has requested rootfinding.


.. c:function:: int ARKodeEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout, sunrealtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param tout: the next time at which a computed solution is desired.
   :param yout: the computed solution vector.
   :param tret: the time corresponding to *yout* (output).
   :param itask: a flag indicating the job of the solver for the next
                 user step.

                 The *ARK_NORMAL* option causes the solver to take internal
                 steps until it has just overtaken a user-specified output
                 time, *tout*, in the direction of integration,
                 i.e. :math:`t_{n-1} <` *tout* :math:`\le t_{n}` for forward
                 integration, or :math:`t_{n} \le` *tout* :math:`< t_{n-1}` for
                 backward integration. If interpolation is enabled (on by
                 default), it will then compute an approximation to the solution
                 :math:`y(tout)` by interpolation (as described in
                 :numref:`ARKODE.Mathematics.Interpolation`). Otherwise, the
                 solution at the time reached by the solver is returned,
                 :math:`y(tret)`.

                 The *ARK_ONE_STEP* option tells the solver to only take a
                 single internal step, :math:`y_{n-1} \to y_{n}`, and return the solution
                 at that point, :math:`y_{n}`, in the vector *yout*.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_ROOT_RETURN: :c:func:`ARKodeEvolve` succeeded, and
                            found one or more roots.  If the number of root functions,
                            *nrtfn*, is greater than 1, call
                            :c:func:`ARKodeGetRootInfo` to see which :math:`g_i` were
                            found to have a root at (*\*tret*).
   :retval ARK_TSTOP_RETURN: :c:func:`ARKodeEvolve` succeeded and
                             returned at *tstop*.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC: ``arkode_mem`` was not allocated.
   :retval ARK_ILL_INPUT: one of the inputs to :c:func:`ARKodeEvolve`
                          is illegal, or some other input to the solver was
                          either illegal or missing.  Details will be
                          provided in the error message.  Typical causes of
                          this failure:

                          (a) A component of the error weight vector became
                              zero during internal time-stepping.

                          (b) The linear solver initialization function (called
                              by the user after calling ``*StepCreate``) failed
                              to set the linear solver-specific *lsolve* field in
                              ``arkode_mem``.

                          (c) A root of one of the root functions was found both at a
                              point :math:`t` and also very near :math:`t`.

                          (d) The initial condition violates the inequality constraints.

   :retval ARK_TOO_MUCH_WORK: the solver took *mxstep* internal steps
                              but could not reach *tout*.  The default value for
                              *mxstep* is *MXSTEP_DEFAULT = 500*.
   :retval ARK_TOO_MUCH_ACC: the solver could not satisfy the accuracy
                             demanded by the user for some internal step.
   :retval ARK_ERR_FAILURE: error test failures occurred either too many
                            times (*ark_maxnef*) during one internal time step
                            or occurred with :math:`|h| = h_{min}`.
   :retval ARK_CONV_FAILURE: either convergence test failures occurred too many
                             times (*ark_maxncf*) during one internal time step
                             or occurred with :math:`|h| = h_{min}`.
   :retval ARK_LINIT_FAIL: the linear solver's initialization function failed.
   :retval ARK_LSETUP_FAIL: the linear solver's setup routine failed in
                            an unrecoverable manner.
   :retval ARK_LSOLVE_FAIL: the linear solver's solve routine failed in
                            an unrecoverable manner.
   :retval ARK_MASSINIT_FAIL: the mass matrix solver's
                              initialization function failed.
   :retval ARK_MASSSETUP_FAIL: the mass matrix solver's setup routine failed.
   :retval ARK_MASSSOLVE_FAIL: the mass matrix solver's solve routine failed.
   :retval ARK_VECTOROP_ERR: a vector operation error occurred.
   :retval ARK_DOMEIG_FAIL: the dominant eigenvalue function failed. It is either
                            not provided or returns an illegal value.
   :retval ARK_MAX_STAGE_LIMIT_FAIL: stepper failed to achieve stable results. Either
                                     reduce the step size or increase the stage_max_limit

   .. note::

      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      ``*StepCreate``.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable.

      All failure return values are negative and so testing the return argument
      for negative values will trap all :c:func:`ARKodeEvolve` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`ARKodeSetStopTime` before the call to
      :c:func:`ARKodeEvolve` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`ARKodeEvolve`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`ARKodeSetStopTime`).

      On any error return in which one or more internal steps were taken
      by :c:func:`ARKodeEvolve`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.OptionalInputs:

Optional input functions
-------------------------

There are numerous optional input parameters that control the behavior
of ARKODE, each of which may be modified from its default value through
calling an appropriate input function.  The following tables list all
optional input functions, grouped by which aspect of ARKODE they control.
Detailed information on the calling syntax and arguments for each
function are then provided following each table.

The optional inputs are grouped into the following categories:

* General ARKODE options (:ref:`ARKODE.Usage.ARKodeInputTable`),
* Step adaptivity solver options (:ref:`ARKODE.Usage.ARKodeAdaptivityInputTable`),
* Implicit stage solver options (:ref:`ARKODE.Usage.ARKodeSolverInputTable`),
* Linear solver interface options (:ref:`ARKODE.Usage.ARKLsInputs`), and
* Rootfinding options (:ref:`ARKODE.Usage.ARKodeRootfindingInputTable`).

For the most casual use of ARKODE, relying on the default set of
solver parameters, the reader can skip to section on user-supplied
functions, :numref:`ARKODE.Usage.UserSupplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so a test on the return arguments for negative values will catch all
errors. Finally, a call to an ``ARKodeSet***`` function can generally be made
from the user's calling program at any time *after* creation of the ARKODE
solver via ``*StepCreate``, and, the function exited successfully, takes effect immediately.
``ARKodeSet***`` functions that cannot be called at any time note
this in the "notes" section of the function documentation.



.. _ARKODE.Usage.ARKodeInputTable:

Optional inputs for ARKODE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=================================================  ==========================================  =======================
Optional input                                     Function name                               Default
=================================================  ==========================================  =======================
Set ARKODE options from the command line or file   :c:func:`ARKodeSetOptions`                  internal
Return ARKODE parameters to their defaults         :c:func:`ARKodeSetDefaults`                 internal
Set integrator method order                        :c:func:`ARKodeSetOrder`                    4
Set dense output interpolation type                :c:func:`ARKodeSetInterpolantType`          stepper-specific
Set dense output polynomial degree                 :c:func:`ARKodeSetInterpolantDegree`        method-dependent
Disable time step adaptivity (fixed-step mode)     :c:func:`ARKodeSetFixedStep`                disabled
Set forward or backward integration direction      :c:func:`ARKodeSetStepDirection`            0.0
Supply an initial step size to attempt             :c:func:`ARKodeSetInitStep`                 estimated
Maximum no. of warnings for :math:`t_n+h = t_n`    :c:func:`ARKodeSetMaxHnilWarns`             10
Maximum no. of internal steps before *tout*        :c:func:`ARKodeSetMaxNumSteps`              500
Maximum absolute step size                         :c:func:`ARKodeSetMaxStep`                  :math:`\infty`
Minimum absolute step size                         :c:func:`ARKodeSetMinStep`                  0.0
Set a value for :math:`t_{stop}`                   :c:func:`ARKodeSetStopTime`                 undefined
Interpolate at :math:`t_{stop}`                    :c:func:`ARKodeSetInterpolateStopTime`      ``SUNFALSE``
Disable the stop time                              :c:func:`ARKodeClearStopTime`               N/A
Supply a pointer for user data                     :c:func:`ARKodeSetUserData`                 ``NULL``
Maximum no. of ARKODE error test failures          :c:func:`ARKodeSetMaxErrTestFails`          7
Set inequality constraints on solution             :c:func:`ARKodeSetConstraints`              ``NULL``
Set max number of constraint failures              :c:func:`ARKodeSetMaxNumConstrFails`        10
Set the checkpointing scheme to use (for adjoint)  :c:func:`ARKodeSetAdjointCheckpointScheme`  ``NULL``
Set the checkpointing step index (for adjoint)     :c:func:`ARKodeSetAdjointCheckpointIndex`   0
Use compensated summation for accumulating time    :c:func:`ARKodeSetUseCompensatedSums`       ``SUNFALSE``
=================================================  ==========================================  =======================


.. c:function:: int ARKodeSetOptions(void* arkode_mem, const char* arkid, const char* file_name, int argc, char* argv[])

   Sets ARKODE options from an array of strings or a file.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param arkid: the prefix for options to read. The default is "arkode".
   :param file_name: the name of a file containing options to read. If this is
                     ``NULL`` or an empty string, ``""``, then no file is read.
   :param argc: length of the ``argv`` array.
   :param argv: an array of strings containing the options to set and their values.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval other: error return value from relevant ARKODE "set" routine.

   **Example usage:**

   In a C or C++ program, the following will enable command-line processing:

   .. code-block:: C

      /* Create ARKODE memory block */
      void* arkode_mem = ARKStepCreate(fe, fi, T0, y, ctx);

      /* Configure ARKODE as normal */
      ...

      /* Override settings with command-line options using default "arkode" prefix */
      flag = ARKodeSetOptions(arkode_mem, NULL, NULL, argc, argv);

   Then when running the program, the user can specify desired options, e.g.,

   .. code-block:: console

      $ ./a.out arkode.order 3 arkode.interpolant_type ARK_INTERP_LAGRANGE

   .. note::

      The ``argc`` and ``argv`` arguments are typically those supplied to the user's
      ``main`` routine however, this is not required. The inputs are left unchanged by
      :c:func:`ARKodeSetOptions`.

      If the ``arkid`` argument is ``NULL``, then the default prefix, ``arkode``, must
      be used for all ARKODE options. Whether ``arkid`` is supplied or not, a ``"."``
      must be used to separate an option key from the prefix.  For example, when
      using the default ``arkid``, the option ``arkode.order`` followed by the value
      can be used to set the method order of accuracy.

      When using a combination of ARKODE integrators (e.g., via MRIStep, SplittingStep,
      or ForcingStep), it is recommended that users call
      :c:func:`ARKodeSetOptions` for each ARKODE integrator using a distinct
      ``arkid`` so they can be controlled separately. For example, "fast" and "slow"
      option prefixes can be used to differentiate between options for the slow and
      fast integrators in an MRI method (i.e., ``fast.order`` and ``slow.order``
      followed by the desired values to set the method order for the fast and slow time
      scales, respectively).

      ARKODE options set via :c:func:`ARKodeSetOptions` will overwrite any previously
      set values. Options are set in the order they are given in ``argv`` and, if an
      option with the same prefix appears multiple times in ``argv``, the value of the
      last occurrence will used.

      The supported option names are noted within the documentation for the
      corresponding ARKODE "set" function. For options that take a
      :c:type:`sunbooleantype` as input, use ``1`` to indicate ``true`` and
      ``0`` for ``false``.

   .. warning::

      This function is not available in the Fortran interface.

      File-based options are not yet supported, so the ``file_name`` argument
      should be set to either ``NULL`` or the empty string ``""``.

   .. versionadded:: 6.5.0


.. c:function:: int ARKodeSetDefaults(void* arkode_mem)

   Resets all optional input parameters to ARKODE's original
   default values.

   :param arkode_mem: pointer to the ARKODE memory block.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      Does not change the *user_data* pointer or any
      parameters within the specified time-stepping module.

      Also leaves alone any data structures or options related to
      root-finding (those can be reset using :c:func:`ARKodeRootInit`).

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the IVP integration method.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ord: requested order of accuracy.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: this option is not supported by the time-stepping module.

   .. note::

      For explicit methods, the allowed values are :math:`2 \le`
      *ord* :math:`\le 8`.  For implicit methods, the allowed values are
      :math:`2\le` *ord* :math:`\le 5`, and for ImEx methods the allowed
      values are :math:`2 \le` *ord* :math:`\le 5`.  Any illegal input
      will result in the default value of 4.

      Since *ord* affects the memory requirements for the internal
      ARKODE memory block, it cannot be changed after the first call to
      :c:func:`ARKodeEvolve`, unless ``*StepReInit`` is called.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.order".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetInterpolantType(void* arkode_mem, int itype)

   Specifies the interpolation type used for dense output (interpolation of
   solution output values) and implicit method predictors. By default,
   Hermite interpolation is used except with SPRK methods where Lagrange
   interpolation is the default.

   This routine must be called *after* the calling a stepper constructor. After
   the first call to :c:func:`ARKodeEvolve` the interpolation type may not be
   changed without first calling a stepper ``ReInit`` function.

   The Hermite interpolation module (``ARK_INTERP_HERMITE``) is described in
   :numref:`ARKODE.Mathematics.Interpolation.Hermite`, and the Lagrange
   interpolation module (``ARK_INTERP_LAGRANGE``) is described in
   :numref:`ARKODE.Mathematics.Interpolation.Lagrange`. ``ARK_INTERP_NONE`` will
   disable interpolation.

   When interpolation is disabled, using rootfinding is not supported, implicit
   methods must use the trivial predictor (the default option), and
   interpolation at stop times cannot be used (interpolating at stop times is
   disabled by default). With interpolation disabled, calling
   :c:func:`ARKodeEvolve` in ``ARK_NORMAL`` mode will return at or past the
   requested output time (setting a stop time may still be used to halt the
   integrator at a specific time).

   Disabling interpolation will reduce the memory footprint of an integrator by
   two or more state vectors (depending on the interpolant type and degree)
   which can be beneficial when interpolation is not needed e.g., when
   integrating to a final time without output in between or using a solver from
   ARKODE as a fast time scale integrator with MRI methods.

   This routine frees any previously-allocated interpolation module, and
   re-creates one according to the specified argument.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param itype: requested interpolant type: ``ARK_INTERP_HERMITE``,
                 ``ARK_INTERP_LAGRANGE``, or ``ARK_INTERP_NONE``

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_MEM_FAIL: the interpolation module could not be allocated.
   :retval ARK_ILL_INPUT: the *itype* argument is not recognized or the
                          interpolation module has already been initialized.

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.interpolant_type".

   .. versionchanged:: 6.1.0

      This function replaces stepper specific versions in ARKStep, ERKStep,
      MRIStep, and SPRKStep.

      Added the ``ARK_INTERP_NONE`` option to disable interpolation.

      Values set by a previous call to :c:func:`ARKStepSetInterpolantDegree` are
      no longer nullified by a call to :c:func:`ARKStepSetInterpolantType`.


.. c:function:: int ARKodeSetInterpolantDegree(void* arkode_mem, int degree)

   Specifies the degree of the polynomial interpolant
   used for dense output (i.e. interpolation of solution output values
   and implicit method predictors).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param degree: requested polynomial degree.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` or the interpolation module are ``NULL``.
   :retval ARK_INTERP_FAIL: this was called after :c:func:`ARKodeEvolve`.
   :retval ARK_ILL_INPUT: an argument had an illegal value or the
                          interpolation module has already been initialized.

   .. note::

      Allowed values are between 0 and 5.

      This routine should be called *before* :c:func:`ARKodeEvolve`. After the
      first call to :c:func:`ARKodeEvolve` the interpolation degree may not be
      changed without first calling ``*StepReInit``.

      If a user calls both this routine and :c:func:`ARKodeSetInterpolantType`, then
      :c:func:`ARKodeSetInterpolantType` must be called first.

      Since the accuracy of any polynomial interpolant is limited by the
      accuracy of the time-step solutions on which it is based, the *actual*
      polynomial degree that is used by ARKODE will be the minimum of
      :math:`q-1` and the input *degree*, for :math:`q > 1` where :math:`q` is
      the order of accuracy for the time integration method.

      When :math:`q=1`, a linear interpolant is the default to ensure values
      obtained by the integrator are returned at the ends of the time
      interval.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.interpolant_degree".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetFixedStep(void* arkode_mem, sunrealtype hfixed)

   Disables time step adaptivity within ARKODE, and specifies the
   fixed time step size to use for the following internal step(s).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hfixed: value of the fixed step size to use.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      Pass 0.0 to return ARKODE to the default (adaptive-step) mode -- this is only
      allowed when using a time-stepping module that supports temporal adaptivity.

      Use of this function is not generally recommended, since it gives no
      assurance of the validity of the computed solutions.  It is
      primarily provided for code-to-code verification testing purposes.

      When using :c:func:`ARKodeSetFixedStep`, any values provided to
      the functions
      :c:func:`ARKodeSetInitStep`,
      :c:func:`ARKodeSetMaxErrTestFails`,
      :c:func:`ARKodeSetCFLFraction`,
      :c:func:`ARKodeSetErrorBias`,
      :c:func:`ARKodeSetFixedStepBounds`,
      :c:func:`ARKodeSetMaxCFailGrowth`,
      :c:func:`ARKodeSetMaxEFailGrowth`,
      :c:func:`ARKodeSetMaxFirstGrowth`,
      :c:func:`ARKodeSetMaxGrowth`,
      :c:func:`ARKodeSetMinReduction`,
      :c:func:`ARKodeSetSafetyFactor`,
      :c:func:`ARKodeSetSmallNumEFails`,
      :c:func:`ARKodeSetStabilityFn`,
      :c:func:`ARKodeSetAdaptController`, and
      :c:func:`ARKodeSetAdaptControllerByName`
      will be ignored, since temporal adaptivity is disabled.

      If both :c:func:`ARKodeSetFixedStep` and
      :c:func:`ARKodeSetStopTime` are used, then the fixed step size
      will be used for all steps until the final step preceding the
      provided stop time (which may be shorter).  To resume use of the
      previous fixed step size, another call to
      :c:func:`ARKodeSetFixedStep` must be made prior to calling
      :c:func:`ARKodeEvolve` to resume integration.

      It is *not* recommended that :c:func:`ARKodeSetFixedStep` be used
      in concert with :c:func:`ARKodeSetMaxStep` or
      :c:func:`ARKodeSetMinStep`, since at best those latter two
      routines will provide no useful information to the solver, and at
      worst they may interfere with the desired fixed step size.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.fixed_step".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetStepDirection(void* arkode_mem, sunrealtype stepdir)

   Specifies the direction of integration (forward or backward).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param stepdir: value whose sign determines the direction. A positive value
                   selects forward integration, a negative value selects
                   backward integration, and zero leaves the current direction
                   unchanged.


   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      The step direction can only be set after a call to either ``*Create``,
      ``*StepReInit``, or :c:func:`ARKodeReset` but before a call to
      :c:func:`ARKodeEvolve`.

      When the direction changes for an adaptive method, the adaptivity
      controller and next step size are reset. A new initial step size will be
      estimated at the next call to :c:func:`ARKodeEvolve` or can be specified
      with :c:func:`ARKodeSetInitStep`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.step_direction".

   .. versionadded:: 6.2.0



.. c:function:: int ARKodeSetInitStep(void* arkode_mem, sunrealtype hin)

   Specifies the initial time step size ARKODE should use after
   initialization, re-initialization, or resetting.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hin: value of the initial step to be attempted :math:`(\ne 0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      Pass 0.0 to use the default value -- this is only
      allowed when using a time-stepping module that supports temporal adaptivity.

      By default, ARKODE estimates the initial step size to be
      :math:`h = \sqrt{\dfrac{2}{\left\| \ddot{y}\right\|}}`, where
      :math:`\ddot{y}` is estimate of the second derivative of the solution
      at :math:`t_0`.

      This routine will also reset the step size and error history.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.init_step".

   .. versionadded:: 6.1.0



.. c:function:: int ARKodeSetMaxHnilWarns(void* arkode_mem, int mxhnil)

   Specifies the maximum number of messages issued by the
   solver to warn that :math:`t+h=t` on the next internal step, before
   ARKODE will instead return with an error.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mxhnil: maximum allowed number of warning messages :math:`(>0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      The default value is 10; set *mxhnil* to zero to specify
      this default.

      A negative value indicates that no warning messages should be issued.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_hnil_warns".

   .. versionadded:: 6.1.0



.. c:function:: int ARKodeSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before ARKODE
   will return with an error.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mxsteps: maximum allowed number of internal steps.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      Passing *mxsteps* = 0 results in ARKODE using the
      default value (500).

      Passing *mxsteps* < 0 disables the test (not recommended).

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_num_steps".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxStep(void* arkode_mem, sunrealtype hmax)

   Specifies the upper bound on the magnitude of the time step size.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hmax: maximum absolute value of the time step size :math:`(\ge 0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Pass *hmax* :math:`\le 0.0` to set the default value of :math:`\infty`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_step".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMinStep(void* arkode_mem, sunrealtype hmin)

   Specifies the lower bound on the magnitude of the time step size.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hmin: minimum absolute value of the time step size :math:`(\ge 0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Pass *hmin* :math:`\le 0.0` to set the default value of 0.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.min_step".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetStopTime(void* arkode_mem, sunrealtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param tstop: stopping time for the integrator.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      The default is that no stop time is imposed.

      Once the integrator returns at a stop time, any future testing for
      ``tstop`` is disabled (and can be re-enabled only though a new call to
      :c:func:`ARKodeSetStopTime`).

      A stop time not reached before a call to ``*StepReInit`` or
      :c:func:`ARKodeReset` will remain active but can be disabled by calling
      :c:func:`ARKodeClearStopTime`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.stop_time".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)

   Specifies that the output solution should be interpolated when the current
   :math:`t` equals the specified ``tstop`` (instead of merely copying the
   internal solution :math:`y_n`).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param interp: flag indicating to use interpolation (1) or copy (0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.interpolate_stop_time".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeClearStopTime(void* arkode_mem)

   Disables the stop time set with :c:func:`ARKodeSetStopTime`.

   :param arkode_mem: pointer to the ARKODE memory block.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      The stop time can be re-enabled though a new call to
      :c:func:`ARKodeSetStopTime`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.clear_stop_time".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* and
   attaches it to the main ARKODE memory block.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param user_data: pointer to the user data.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      If specified, the pointer to *user_data* is passed to all
      user-supplied functions for which it is an argument; otherwise
      ``NULL`` is passed.

      If *user_data* is needed in user preconditioner functions, the call to
      this function must be made *before* any calls to
      :c:func:`ARKodeSetLinearSolver` and/or :c:func:`ARKodeSetMassLinearSolver`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxErrTestFails(void* arkode_mem, int maxnef)

   Specifies the maximum number of error test failures
   permitted in attempting one step, before returning with an error.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param maxnef: maximum allowed number of error test failures :math:`(>0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      The default value is 7; set *maxnef* :math:`\le 0`
      to specify this default.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_err_test_fails".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetConstraints(void* arkode_mem, N_Vector constraints)

   Specifies a vector defining inequality constraints for each component of the
   solution vector :math:`y`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param constraints: vector of constraint flags. Each component specifies
                       the type of solution constraint:

                       .. math::

                          \texttt{constraints[i]} = \left\{ \begin{array}{rcl}
                          0.0  &\Rightarrow\;& \text{no constraint is imposed on}\; y_i,\\
                          1.0  &\Rightarrow\;& y_i \geq 0,\\
                          -1.0  &\Rightarrow\;& y_i \leq 0,\\
                          2.0  &\Rightarrow\;& y_i > 0,\\
                          -2.0  &\Rightarrow\;& y_i < 0.\\
                          \end{array}\right.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: the constraints vector contains illegal values.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      The presence of a non-``NULL`` constraints vector that is not 0.0
      in all components will cause constraint checking to be performed. However, a
      call with 0.0 in all components of ``constraints`` will result in an illegal
      input return. A ``NULL`` constraints vector will disable constraint checking.

      After a call to :c:func:`ARKodeResize` inequality constraint checking
      will be disabled and a call to :c:func:`ARKodeSetConstraints` is
      required to re-enable constraint checking.

      Since constraint-handling is performed through cutting time steps that would
      violate the constraints, it is possible that this feature will cause some
      problems to fail due to an inability to enforce constraints even at the
      minimum time step size.  Additionally, the features :c:func:`ARKodeSetConstraints`
      and :c:func:`ARKodeSetFixedStep` are incompatible, and should not be used
      simultaneously.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxNumConstrFails(void* arkode_mem, int maxfails)

   Specifies the maximum number of constraint failures in a step before ARKODE
   will return with an error.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param maxfails: maximum allowed number of constrain failures.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Passing *maxfails* <= 0 results in ARKODE using the
      default value (10).

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_num_constr_fails".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetAdjointCheckpointScheme(void* arkode_mem, SUNAdjointCheckpointScheme checkpoint_scheme)

   Specifies the :c:type:`SUNAdjointCheckpointScheme` to use for saving states
   during the forward integration, and loading states during backward integration
   of an adjoint system.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param checkpoint_scheme: the checkpoint scheme to use, or ``NULL`` to disable checkpointing.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.3.0


.. c:function:: int ARKodeSetAdjointCheckpointIndex(void* arkode_mem, suncountertype step_index)

   Specifies the step index (that is step number) to insert the next checkpoint at.

   This is incremented along with the step count, but it is useful to be able to reset
   this index during recomputations of missing states during the backward adjoint integration.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param step_idx: the step to insert the next checkpoint at.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.3.0


.. c:function:: int ARKodeSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)

   Specifies if compensated summations should be used within ARKODE where supported.
   Currently, all ARKODE modules support compensated summation for accumulating time.

   SPRKStep also supports an alternative stepping algorithm based on compensated
   summation which will be enabled/disabled by this function. This increases the
   computational cost by 2 extra vector operations per stage and an additional
   5 per time step. It also requires one extra vector to be stored. However, it
   is significantly more robust to roundoff error accumulation.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param onoff: should compensated summation be used (1) or not (0)

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the ARKODE memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.use_compensated_sums".




.. _ARKODE.Usage.ARKodeAdaptivityInputTable:

Optional inputs for time step adaptivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of ARKODE's time step adaptivity
algorithm, including how each of the parameters below is used within
the code, is provided in :numref:`ARKODE.Mathematics.Adaptivity`.


.. cssclass:: table-bordered

=========================================================   ==========================================  ========
Optional input                                              Function name                               Default
=========================================================   ==========================================  ========
Provide a :c:type:`SUNAdaptController` for ARKODE to use    :c:func:`ARKodeSetAdaptController`          I
Specify a :c:type:`SUNAdaptController` for ARKODE to use    :c:func:`ARKodeSetAdaptControllerByName`    I
Adjust the method order used in the controller              :c:func:`ARKodeSetAdaptivityAdjustment`     0
Explicit stability safety factor                            :c:func:`ARKodeSetCFLFraction`              0.5
Time step error bias factor                                 :c:func:`ARKodeSetErrorBias`                1.0
Bounds determining no change in step size                   :c:func:`ARKodeSetFixedStepBounds`          1.0  1.0
Maximum step growth factor on convergence fail              :c:func:`ARKodeSetMaxCFailGrowth`           0.25
Maximum step growth factor on error test fail               :c:func:`ARKodeSetMaxEFailGrowth`           0.3
Maximum first step growth factor                            :c:func:`ARKodeSetMaxFirstGrowth`           10000.0
Maximum allowed general step growth factor                  :c:func:`ARKodeSetMaxGrowth`                20.0
Minimum allowed step reduction factor on error test fail    :c:func:`ARKodeSetMinReduction`             0.1
Time step safety factor                                     :c:func:`ARKodeSetSafetyFactor`             0.9
Error fails before ``MaxEFailGrowth`` takes effect          :c:func:`ARKodeSetSmallNumEFails`           2
Explicit stability function                                 :c:func:`ARKodeSetStabilityFn`              none
Set accumulated error estimation type                       :c:func:`ARKodeSetAccumulatedErrorType`     none
Reset accumulated error                                     :c:func:`ARKodeResetAccumulatedError`
=========================================================   ==========================================  ========



.. c:function:: int ARKodeSetAdaptController(void* arkode_mem, SUNAdaptController C)

   Sets a user-supplied time-step controller object.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param C: user-supplied time adaptivity controller.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_MEM_FAIL: *C* was ``NULL`` and the I controller could not be allocated.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      If *C* is ``NULL`` then the I controller will be created (see :numref:`SUNAdaptController.Soderlind`).

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Not all time-stepping modules are compatible with all types of :c:type:`SUNAdaptController`
      objects.  While all steppers that support temporal adaptivity support controllers with
      :c:type:`SUNAdaptController_Type` type ``SUN_ADAPTCONTROLLER_H``, only MRIStep supports
      inputs with type ``SUN_ADAPTCONTROLLER_MRI_H_TOL``.

   .. versionadded:: 6.1.0

   .. versionchanged:: 6.3.0

      The default controller was changed from PID to I.




.. c:function:: int ARKodeSetAdaptControllerByName(void* arkode_mem, const char* cname)

   Sets a user-supplied time step controller object by name.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param cname: name of the time adaptivity controller to use.  Allowable values
                 currently include ``"Soderlind"``, ``"PID"``, ``"PI"``, ``"I"``,
                 ``"ExpGus"``, ``"ImpGus"``, ``"ImExGus"``, ``"H0211"``, ``"H0321"``,
                 ``"H211"``, and ``"H312"``. For information on these options, see
                 :numref:`SUNAdaptController.Soderlind` and
                 :numref:`SUNAdaptController.ImExGus`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_ILL_INPUT: ``cname`` did not match an allowed value.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      It is not possible to adjust the internal controller parameters when using this
      function.  Users who wish to adjust these parameters should create and configure
      the :c:type:`SUNAdaptController` object manually, and then call
      :c:func:`ARKodeSetAdaptController`.

   .. versionadded:: 6.3.0


.. c:function:: int ARKodeSetAdaptivityAdjustment(void* arkode_mem, int adjust)

   Called by a user to adjust the method order supplied to the temporal adaptivity
   controller.  For example, if the user expects order reduction due to problem stiffness,
   they may request that the controller assume a reduced order of accuracy for the method
   by specifying a value :math:`adjust < 0`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param adjust: adjustment factor (default is 0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      This should be called prior to calling :c:func:`ARKodeEvolve`, and can only be
      reset following a call to ``*StepReInit``.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.adaptivity_adjustment".

   .. versionadded:: 6.1.0

   .. versionchanged:: 6.3.0

      The default value was changed from -1 to 0


.. c:function:: int ARKodeSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac)

   Specifies the fraction of the estimated explicitly stable step to use.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param cfl_frac: maximum allowed fraction of explicitly stable step (default is 0.5).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any non-positive parameter will imply a reset to the default
      value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.cfl_fraction".

   .. versionadded:: 6.1.0

   .. versionchanged:: 6.4.0

      The restriction that ``cfl_frac`` is less than one has been removed.


.. c:function:: int ARKodeSetErrorBias(void* arkode_mem, sunrealtype bias)

   Specifies the bias to be applied to the error estimates within
   accuracy-based adaptivity strategies.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param bias: bias applied to error in accuracy-based time
                step estimation (default is 1.0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value below 1.0 will imply a reset to the default value.

      If both this and one of the stepper ``SetAdaptivityMethod`` functions or
      :c:func:`ARKodeSetAdaptController` will be called, then this routine must be called
      *second*.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.error_bias".

   .. versionadded:: 6.1.0

   .. versionchanged:: 6.3.0

      The default value was changed from 1.5 to 1.0


.. c:function:: int ARKodeSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub)

   Specifies the step growth interval in which the step size will remain unchanged.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lb: lower bound on window to leave step size fixed (default is 1.0).
   :param ub: upper bound on window to leave step size fixed (default is 1.0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any interval *not* containing 1.0 will imply a reset to the default values.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.fixed_step_bounds".

   .. versionadded:: 6.1.0

   .. versionchanged:: 6.3.0

      The default upper bound was changed from 1.5 to 1.0


.. c:function:: int ARKodeSetMaxCFailGrowth(void* arkode_mem, sunrealtype etacf)

   Specifies the maximum step size growth factor upon an algebraic
   solver convergence failure on a stage solve within a step, :math:`\eta_{cf}` from
   :numref:`ARKODE.Mathematics.Error.Nonlinear`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param etacf: time step reduction factor on a nonlinear solver
                 convergence failure (default is 0.25).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value outside the interval :math:`(0,1]` will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_cfail_growth".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf)

   Specifies the maximum step size growth factor upon multiple successive
   accuracy-based error failures in the solver.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param etamxf: time step reduction factor on multiple error fails (default is 0.3).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value outside the interval :math:`(0,1]` will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_efail_growth".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1)

   Specifies the maximum allowed growth factor in step size following the very
   first integration step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param etamx1: maximum allowed growth factor after the first time
                  step (default is 10000.0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value :math:`\le 1.0` will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_first_growth".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth)

   Specifies the maximum allowed growth factor in step size between
   consecutive steps in the integration process.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mx_growth: maximum allowed growth factor between consecutive time steps (default is 20.0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value :math:`\le 1.0` will imply a reset to the default
      value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_growth".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMinReduction(void* arkode_mem, sunrealtype eta_min)

   Specifies the minimum allowed reduction factor in step size between
   step attempts, resulting from a temporal error failure in the integration
   process.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param eta_min: minimum allowed reduction factor in time step after an error
                   test failure (default is 0.1).
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value outside the interval :math:`(0,1)` will imply a reset to
      the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.min_reduction".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetSafetyFactor(void* arkode_mem, sunrealtype safety)

   Specifies the safety factor to be applied to the accuracy-based
   estimated step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param safety: safety factor applied to accuracy-based time step (default is 0.9).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value :math:`\le 0` will imply a reset to the default
      value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.safety_factor".

   .. versionadded:: 6.1.0

   .. versionchanged:: 6.3.0

      The default default was changed from 0.96 to 0.9. The maximum value is now
      exactly 1.0 rather than strictly less than 1.0.


.. c:function:: int ARKodeSetSmallNumEFails(void* arkode_mem, int small_nef)

   Specifies the threshold for "multiple" successive error failures
   before the *etamxf* parameter from
   :c:func:`ARKodeSetMaxEFailGrowth` is applied.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param small_nef: bound to determine 'multiple' for *etamxf* (default is 2).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      Any value :math:`\le 0` will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.small_num_efails".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)

   Sets the problem-dependent function to estimate a stable
   time step size for the explicit portion of the ODE system.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param EStab: name of user-supplied stability function.
   :param estab_data: pointer to user data passed to *EStab* every time
                      it is called.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

      This function should return an estimate of the absolute
      value of the maximum stable time step for the explicit portion of
      the ODE system.  It is not required, since accuracy-based
      adaptivity may be sufficient for retaining stability, but this can
      be quite useful for problems where the explicit right-hand side
      function :math:`f^E(t,y)` contains stiff terms.

   .. versionadded:: 6.1.0


The following routines are used to control algorithms that ARKODE can use to estimate
the accumulated temporal error over multiple time steps.  While these may be informational
for users on their applications, this functionality is required when using multirate
temporal adaptivity in MRIStep via the :ref:`SUNAdaptController_MRIHTol <SUNAdaptController.MRIHTol>`
module.  For time-stepping modules that compute both a solution and embedding, :math:`y_n`
and :math:`\tilde{y}_n`, these may be combined to create a vector-valued local temporal error
estimate for the current internal step, :math:`y_n - \tilde{y}_n`.  These local errors may be
accumulated by ARKODE in a variety of ways, as determined by the enumerated type
:c:enum:`ARKAccumError`.  In each of the cases below, the accumulation is taken over all steps
since the most recent call to either :c:func:`ARKodeSetAccumulatedErrorType` or
:c:func:`ARKodeResetAccumulatedError`. Below the set :math:`\mathcal{S}` contains
the indices of the steps since the last call to either of the aforementioned functions.
The norm is taken using the tolerance-informed error-weight vector (see
:c:func:`ARKodeGetErrWeights`), and ``reltol`` is the user-specified relative solution
tolerance.

.. c:enum:: ARKAccumError

   The type of error accumulation that ARKODE should use.

   .. versionadded:: 6.2.0

   .. c:enumerator:: ARK_ACCUMERROR_NONE

      No accumulation should be performed

   .. c:enumerator:: ARK_ACCUMERROR_MAX

      Computes :math:`\text{reltol} \max\limits_{i \in \mathcal{S}} \|y_i - \tilde{y}_i\|_{WRMS}`

   .. c:enumerator:: ARK_ACCUMERROR_SUM

      Computes :math:`\text{reltol} \sum\limits_{i \in \mathcal{S}} \|y_i - \tilde{y}_i\|_{WRMS}`

   .. c:enumerator:: ARK_ACCUMERROR_AVG

      Computes :math:`\frac{\text{reltol}}{\Delta t_{\mathcal{S}}} \sum\limits_{i \in \mathcal{S}} h_i \|y_i - \tilde{y}_i\|_{WRMS}`,
      where :math:`h_i` is the step size used when computing :math:`y_i`, and
      :math:`\Delta t_{\mathcal{S}}` denotes the elapsed time over which
      :math:`\mathcal{S}` is taken.


.. c:function:: int ARKodeSetAccumulatedErrorType(void* arkode_mem, ARKAccumError accum_type)

   Sets the strategy to use for accumulating a temporal error estimate
   over multiple time steps.  By default, ARKODE will not accumulate any
   local error estimates (i.e., the default *accum_type* is ``ARK_ACCUMERROR_NONE``).

   A non-default error accumulation strategy can be disabled by calling
   :c:func:`ARKodeSetAccumulatedErrorType` with the argument ``ARK_ACCUMERROR_NONE``.

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.accumulated_error_type".

   :param arkode_mem: pointer to the ARKODE memory block.
   :param accum_type: accumulation strategy.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_STEPPER_UNSUPPORTED: temporal error estimation is not supported
                                    by the current time-stepping module.

   .. versionadded:: 6.2.0


.. c:function:: int ARKodeResetAccumulatedError(void* arkode_mem)

   Resets the accumulated temporal error estimate, that was triggered by a previous call to
   :c:func:`ARKodeSetAccumulatedErrorType`.

   :param arkode_mem: pointer to the ARKODE memory block.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_STEPPER_UNSUPPORTED: temporal error estimation is not supported
                                    by the current time-stepping module.

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.reset_accumulated_error".

   .. versionadded:: 6.2.0



.. _ARKODE.Usage.ARKodeSolverInputTable:

Optional inputs for implicit stage solves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation for the nonlinear solver strategies used
by ARKODE, including how each of the parameters below is used within
the code, is provided in :numref:`ARKODE.Mathematics.Nonlinear`.


.. cssclass:: table-bordered

==============================================================  ======================================  ============
Optional input                                                  Function name                           Default
==============================================================  ======================================  ============
Specify that the implicit RHS is linear                         :c:func:`ARKodeSetLinear`               ``SUNFALSE``
Specify that the implicit RHS nonlinear                         :c:func:`ARKodeSetNonlinear`            ``SUNTRUE``
Specify that the implicit RHS is autonomous                     :c:func:`ARKodeSetAutonomous`           ``SUNFALSE``
Implicit predictor method                                       :c:func:`ARKodeSetPredictorMethod`      0
User-provided implicit stage predictor                          :c:func:`ARKodeSetStagePredictFn`       ``NULL``
RHS function for nonlinear system evaluations                   :c:func:`ARKodeSetNlsRhsFn`             ``NULL``
Maximum number of nonlinear iterations                          :c:func:`ARKodeSetMaxNonlinIters`       3
Coefficient in the nonlinear convergence test                   :c:func:`ARKodeSetNonlinConvCoef`       0.1
Nonlinear convergence rate constant                             :c:func:`ARKodeSetNonlinCRDown`         0.3
Nonlinear residual divergence ratio                             :c:func:`ARKodeSetNonlinRDiv`           2.3
Maximum number of convergence failures                          :c:func:`ARKodeSetMaxConvFails`         10
Specify if the implicit RHS is deduced after a nonlinear solve  :c:func:`ARKodeSetDeduceImplicitRhs`    ``SUNFALSE``
==============================================================  ======================================  ============





.. c:function:: int ARKodeSetLinear(void* arkode_mem, int timedepend)

   Specifies that the implicit portion of the problem is linear.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param timedepend: flag denoting whether the Jacobian of
                      :math:`f^I(t,y)` is time-dependent (1) or not (0).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Tightens the linear solver tolerances and takes only a
      single Newton iteration.  Calls :c:func:`ARKodeSetDeltaGammaMax`
      to enforce Jacobian recomputation when the step size ratio changes
      by more than 100 times the unit roundoff (since nonlinear
      convergence is not tested).  Only applicable when used in
      combination with the modified or inexact Newton iteration (not the
      fixed-point solver).

      When :math:`f^I(t,y)` is time-dependent, all linear solver structures
      (Jacobian, preconditioner) will be updated preceding *each* implicit
      stage.  Thus one must balance the relative costs of such recomputation
      against the benefits of requiring only a single Newton linear solve.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetNonlinear(void* arkode_mem)

   Specifies that the implicit portion of the problem is nonlinear.

   :param arkode_mem: pointer to the ARKODE memory block.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is the default behavior of ARKODE, so the function
      is primarily useful to undo a previous call to
      :c:func:`ARKodeSetLinear`.  Calls
      :c:func:`ARKodeSetDeltaGammaMax` to reset the step size ratio
      threshold to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.nonlinear".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetAutonomous(void* arkode_mem, sunbooleantype autonomous)

   Specifies that the implicit portion of the problem is autonomous i.e., does
   not explicitly depend on time.

   When using an implicit or ImEx method with the trivial predictor, this option
   enables reusing the implicit right-hand side evaluation at the predicted
   state across stage solves within a step. This reuse reduces the total number
   of implicit RHS function evaluations.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param autonomous: flag denoting if the implicit RHS function,
                      :math:`f^I(t,y)`, is autonomous (``SUNTRUE``) or
                      non-autonomous (``SUNFALSE``).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. warning::

      Results may differ when enabling both :c:func:`ARKodeSetAutonomous` and
      :c:func:`ARKodeSetDeduceImplicitRhs` with a stiffly accurate implicit
      method and using the trivial predictor. The differences are due to reusing
      the deduced implicit right-hand side (RHS) value in the initial nonlinear
      residual computation rather than evaluating the implicit RHS function. The
      significance of the difference will depend on how well the deduced RHS
      approximates the RHS evaluated at the trivial predictor. This behavior can
      be observed in ``examples/arkode/C_serial/ark_brusselator.c`` by comparing
      the outputs with :c:func:`ARKodeSetAutonomous` enabled/disabled.

      Similarly programs that assume the nonlinear residual will always call the
      implicit RHS function will need to be updated to account for the RHS value
      reuse when using :c:func:`ARKodeSetAutonomous`. For example,
      ``examples/arkode/C_serial/ark_KrylovDemo_prec.c`` assumes that the
      nonlinear residual will be called and will evaluate the implicit RHS
      function before calling the preconditioner setup function. Based on this
      assumption, this example code saves some computations in the RHS
      evaluation for reuse in the preconditioner setup. However, when
      :c:func:`ARKodeSetAutonomous` is enabled, the call to the nonlinear
      residual before the preconditioner setup reuses a saved RHS evaluation and
      the saved data is actually from an earlier RHS evaluation that is not
      consistent with the state and RHS values passed to the preconditioner
      setup function. For this example, the code should not save data in the RHS
      evaluation but instead evaluate the necessary quantities within the
      preconditioner setup function using the input values.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.autonomous".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetPredictorMethod(void* arkode_mem, int method)

   Specifies the method from :numref:`ARKODE.Mathematics.Predictors` to use
   for predicting implicit solutions.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param method: method choice (0 :math:`\le` *method* :math:`\le` 4):

                  * 0 is the trivial predictor,

                  * 1 is the maximum order (dense output) predictor,

                  * 2 is the variable order predictor, that decreases the
                    polynomial degree for more distant RK stages,

                  * 3 is the cutoff order predictor, that uses the maximum order
                    for early RK stages, and a first-order predictor for distant
                    RK stages,

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default value is 0.  If *method* is set to an
      undefined value, this default predictor will be used.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.predictor_method".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage)

   Sets the user-supplied function to update the implicit stage predictor prior to
   execution of the nonlinear or linear solver algorithms that compute the implicit stage solution.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param PredictStage: name of user-supplied predictor function.  If ``NULL``, then any
                        previously-provided stage prediction function will be disabled.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      See :numref:`ARKODE.Usage.StagePredictFn` for more information on
      this user-supplied routine.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fi)

   Specifies an alternative implicit right-hand side function for evaluating
   :math:`f^I(t,y)` within nonlinear system function evaluations
   :eq:`ARKODE_Residual_MeqI` - :eq:`ARKODE_Residual_MTimeDep`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nls_fi: the alternative C function for computing the right-hand side
                  function :math:`f^I(t,y)` in the ODE.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default is to use the implicit right-hand side function
      provided to the stepper constructor in nonlinear system functions. If the
      input implicit right-hand side function is ``NULL``, the default is used.

      When using a non-default nonlinear solver, this function must be called
      *after* :c:func:`ARKodeSetNonlinearSolver`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxNonlinIters(void* arkode_mem, int maxcor)

   Specifies the maximum number of nonlinear solver
   iterations permitted per implicit stage solve within each time step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param maxcor: maximum allowed solver iterations per stage :math:`(>0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value or if the SUNNONLINSOL module is ``NULL``.
   :retval ARK_NLS_OP_ERR: the SUNNONLINSOL object returned a failure flag.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default value is 3; set *maxcor* :math:`\le 0`
      to specify this default.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_nonlin_iters".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef)

   Specifies the safety factor :math:`\epsilon` used within the nonlinear
   solver convergence test :eq:`ARKODE_NonlinearTolerance`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nlscoef: coefficient in nonlinear solver convergence test :math:`(>0.0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default value is 0.1; set *nlscoef* :math:`\le 0`
      to specify this default.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.nonlin_conv_coef".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetNonlinCRDown(void* arkode_mem, sunrealtype crdown)

   Specifies the constant :math:`c_r` used in estimating the nonlinear solver convergence rate :eq:`ARKODE_NonlinearCRate`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param crdown: nonlinear convergence rate estimation constant (default is 0.3).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Any non-positive parameter will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.nonlin_crdown".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv)

   Specifies the nonlinear correction threshold :math:`r_{div}` from
   :eq:`ARKODE_NonlinearDivergence`, beyond which the iteration will be declared divergent.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rdiv: tolerance on nonlinear correction size ratio to
                declare divergence (default is 2.3).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Any non-positive parameter will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.nonlin_rdiv".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMaxConvFails(void* arkode_mem, int maxncf)

   Specifies the maximum number of nonlinear solver convergence
   failures permitted during one step, :math:`max_{ncf}` from
   :numref:`ARKODE.Mathematics.Error.Nonlinear`, before ARKODE will return with
   an error.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param maxncf: maximum allowed nonlinear solver convergence failures
                  per step :math:`(>0)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default value is 10; set *maxncf* :math:`\le 0`
      to specify this default.

      Upon each convergence failure, ARKODE will first call the Jacobian
      setup routine and try again (if a Newton method is used).  If a
      convergence failure still occurs, the time step size is reduced by
      the factor *etacf* (set within :c:func:`ARKodeSetMaxCFailGrowth`).

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_conv_fails".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetDeduceImplicitRhs(void *arkode_mem, sunbooleantype deduce)

   Specifies if implicit stage derivatives are deduced without evaluating
   :math:`f^I`. See :numref:`ARKODE.Mathematics.Nonlinear` for more details.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param deduce: if ``SUNFALSE`` (default), the stage derivative is obtained
                  by evaluating :math:`f^I` with the stage solution returned from the
                  nonlinear solver. If ``SUNTRUE``, the stage derivative is deduced
                  without an additional evaluation of :math:`f^I`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.deduce_implicit_rhs".

   .. versionadded:: 6.1.0


.. _ARKODE.Usage.ARKLsInputs:


Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of the linear solver methods
available to ARKODE is provided in :numref:`ARKODE.Mathematics.Linear`.  We group
the user-callable routines into four categories: general routines concerning
the update frequency for matrices and/or preconditioners, optional inputs for
matrix-based linear solvers, optional inputs for matrix-free linear solvers,
and optional inputs for iterative linear solvers.  We note that the
matrix-based and matrix-free groups are mutually exclusive, whereas the
"iterative" tag can apply to either case.



.. _ARKODE.Usage.ARKLsInputs.General:

.. index::
   single: optional input; generic linear solver interface (ARKODE)

Optional inputs for the ARKLS linear solver interface
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As discussed in :numref:`ARKODE.Mathematics.Linear.Setup`, ARKODE
strives to reuse matrix and preconditioner data for as many solves as
possible to amortize the high costs of matrix construction and
factorization.  To that end, ARKODE provides user-callable
routines to modify this behavior.  Recall that the
Newton system matrices that arise within an implicit stage solve are
:math:`\mathcal{A}(t,z) \approx M(t) - \gamma J(t,z)`, where the
implicit right-hand side function has Jacobian matrix
:math:`J(t,z) = \frac{\partial f^I(t,z)}{\partial z}`.

The matrix or preconditioner for :math:`\mathcal{A}` can only be
updated within a call to the linear solver "setup" routine.  In
general, the frequency with which the linear solver setup routine is
called may be controlled with the *msbp* argument to
:c:func:`ARKodeSetLSetupFrequency`.  When this occurs, the
validity of :math:`\mathcal{A}` for successive time steps
intimately depends on whether the corresponding :math:`\gamma` and
:math:`J` inputs remain valid.

At each call to the linear solver setup routine the decision to update
:math:`\mathcal{A}` with a new value of :math:`\gamma`, and to reuse
or reevaluate Jacobian information, depends on several factors including:

* the success or failure of previous solve attempts,
* the success or failure of the previous time step attempts,
* the change in :math:`\gamma` from the value used when constructing :math:`\mathcal{A}`, and
* the number of steps since Jacobian information was last evaluated.

Jacobian information is considered out-of-date when :math:`msbj` or more steps
have been completed since the last update, in which case it will be recomputed during the next
linear solver setup call. The value of :math:`msbj` is controlled with the
``msbj`` argument to :c:func:`ARKodeSetJacEvalFrequency`.

For linear-solvers with user-supplied preconditioning the above factors are used
to determine whether to recommend updating the Jacobian information in the
preconditioner (i.e., whether to set *jok* to ``SUNFALSE`` in calling the
user-supplied :c:type:`ARKLsPrecSetupFn`). For matrix-based linear solvers
these factors determine whether the matrix :math:`J(t,y) = \frac{\partial f^I(t,y)}{\partial y}`
should be updated (either with an internal finite difference approximation or
a call to the user-supplied :c:type:`ARKLsJacFn`); if not then the previous
value is reused and the system matrix :math:`\mathcal{A}(t,y) \approx M(t) - \gamma J(t,y)`
is recomputed using the current :math:`\gamma` value.



.. _ARKODE.Usage.ARKLsInputs.General.Table:
.. table:: Optional inputs for the ARKLS linear solver interface

   =============================================  ====================================  ============
   Optional input                                 Function name                         Default
   =============================================  ====================================  ============
   Max change in step signaling new :math:`J`     :c:func:`ARKodeSetDeltaGammaMax`      0.2
   Linear solver setup frequency                  :c:func:`ARKodeSetLSetupFrequency`    20
   Jacobian / preconditioner update frequency     :c:func:`ARKodeSetJacEvalFrequency`   51
   =============================================  ====================================  ============


.. c:function:: int ARKodeSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax)

   Specifies a scaled step size ratio tolerance, :math:`\Delta\gamma_{max}` from
   :numref:`ARKODE.Mathematics.Linear.Setup`, beyond which the linear solver
   setup routine will be signaled.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param dgmax: tolerance on step size ratio change before calling
                 linear solver setup routine (default is 0.2).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Any non-positive parameter will imply a reset to the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.delta_gamma_max".

   .. versionadded:: 6.1.0

.. index::
   single: optional input; linear solver setup frequency (ARKODE)

.. c:function:: int ARKodeSetLSetupFrequency(void* arkode_mem, int msbp)

   Specifies the frequency of calls to the linear solver setup
   routine, :math:`msbp` from :numref:`ARKODE.Mathematics.Linear.Setup`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param msbp: the linear solver setup frequency.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Positive values of **msbp** specify the linear solver setup frequency. For
      example, an input of 1 means the setup function will be called every time
      step while an input of 2 means it will be called called every other time
      step. If **msbp** is 0, the default value of 20 will be used. A negative
      value forces a linear solver step at each implicit stage.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.lsetup_frequency".

   .. versionadded:: 6.1.0


.. index::
   single: optional input; Jacobian update frequency (ARKODE)
   single: optional input; preconditioner update frequency (ARKODE)

.. c:function:: int ARKodeSetJacEvalFrequency(void* arkode_mem, long int msbj)

   Specifies the number of steps after which the Jacobian information is
   considered out-of-date, :math:`msbj` from :numref:`ARKODE.Mathematics.Linear.Setup`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param msbj: the Jacobian re-computation or preconditioner update frequency.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      If ``nstlj`` is the step number at which the Jacobian information was
      lasted updated and ``nst`` is the current step number,
      ``nst - nstlj >= msbj`` indicates that the Jacobian information will be updated
      during the next linear solver setup call.

      As the Jacobian update frequency is only checked *within* calls to the
      linear solver setup routine, Jacobian information may be more than
      ``msbj`` steps old when updated depending on when a linear solver setup
      call occurs. See :numref:`ARKODE.Mathematics.Linear.Setup`
      for more information on when linear solver setups are performed.

      Passing a value *msbj* :math:`\le 0` indicates to use the
      default value of 51.

      This function must be called *after* the ARKLS system solver interface has
      been initialized through a call to :c:func:`ARKodeSetLinearSolver`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.jac_eval_frequency".

   .. versionadded:: 6.1.0






.. _ARKODE.Usage.ARKLsInputs.MatrixBased:

Optional inputs for matrix-based ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

=========================================  ========================================  =============
Optional input                             Function name                             Default
=========================================  ========================================  =============
Jacobian function                          :c:func:`ARKodeSetJacFn`                  ``DQ``
Linear system function                     :c:func:`ARKodeSetLinSysFn`               internal
Mass matrix function                       :c:func:`ARKodeSetMassFn`                 none
Enable or disable linear solution scaling  :c:func:`ARKodeSetLinearSolutionScaling`  on
=========================================  ========================================  =============

When using matrix-based linear solver modules, the ARKLS solver interface needs
a function to compute an approximation to the Jacobian matrix :math:`J(t,y)` or
the linear system :math:`\mathcal{A}(t,y) = M(t) - \gamma J(t,y)`.

For :math:`J(t,y)`, the ARKLS interface is packaged with a routine that can approximate
:math:`J` if the user has selected either the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
:ref:`SUNMATRIX_BAND <SUNMatrix.Band>` objects.  Alternatively,
the user can supply a custom Jacobian function of type :c:func:`ARKLsJacFn` -- this is
*required* when the user selects other matrix formats.  To specify a user-supplied
Jacobian function, ARKODE provides the function :c:func:`ARKodeSetJacFn`.

Alternatively, a function of type :c:func:`ARKLsLinSysFn` can be provided to
evaluate the matrix :math:`\mathcal{A}(t,y)`. By default, ARKLS uses an
internal linear system function leveraging the SUNMATRIX API to form the matrix
:math:`\mathcal{A}(t,y)` by combining the matrices :math:`M(t)` and :math:`J(t,y)`.
To specify a user-supplied linear system function instead, ARKODE provides the function
:c:func:`ARKodeSetLinSysFn`.

If the ODE system involves a non-identity mass matrix, :math:`M\ne I`, matrix-based linear
solver modules require a function to compute an approximation to the mass matrix :math:`M(t)`.
There is no default difference quotient approximation (for any matrix type), so this
routine must be supplied by the user. This function must be of type
:c:func:`ARKLsMassFn`, and should be set using the function
:c:func:`ARKodeSetMassFn`.

In either case (:math:`J(t,y)` versus :math:`\mathcal{A}(t,y)` is supplied) the matrix
information will be updated infrequently to reduce matrix construction and, with direct
solvers, factorization costs. As a result the value of :math:`\gamma` may not be current
and a scaling factor is applied to the solution of the linear system to account for
the lagged value of :math:`\gamma`. See :numref:`SUNLinSol.Lagged_matrix` for more details.
The function :c:func:`ARKodeSetLinearSolutionScaling` can be used to disable this
scaling when necessary, e.g., when providing a custom linear solver that updates the
matrix using the current :math:`\gamma` as part of the solve.

The ARKLS interface passes the user data pointer to the Jacobian, linear
system, and mass matrix functions. This allows the user to create an arbitrary
structure with relevant problem data and access it during the execution of the
user-supplied Jacobian, linear system or mass matrix functions, without using global
data in the program. The user data pointer may be specified through
:c:func:`ARKodeSetUserData`.



.. c:function:: int ARKodeSetJacFn(void* arkode_mem, ARKLsJacFn jac)

   Specifies the Jacobian approximation routine to
   be used for the matrix-based solver with the ARKLS interface.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param jac: name of user-supplied Jacobian approximation function.

   :retval ARKLS_SUCCESS:  the function exited successfully.
   :retval ARKLS_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This routine must be called after the ARKLS linear
      solver interface has been initialized through a call to
      :c:func:`ARKodeSetLinearSolver`.

      By default, ARKLS uses an internal difference quotient function for
      the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and
      :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` modules.  If ``NULL`` is passed
      in for *jac*, this default is used. An error will occur if no *jac* is
      supplied when using other matrix types.

      The function type :c:func:`ARKLsJacFn` is described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys)

   Specifies the linear system approximation routine to be used for the
   matrix-based solver with the ARKLS interface.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param linsys: name of user-supplied linear system approximation function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This routine must be called after the ARKLS linear
      solver interface has been initialized through a call to
      :c:func:`ARKodeSetLinearSolver`.

      By default, ARKLS uses an internal linear system function that leverages the
      SUNMATRIX API to form the system :math:`M - \gamma J`.  If ``NULL`` is passed
      in for *linsys*, this default is used.

      The function type :c:func:`ARKLsLinSysFn` is described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMassFn(void* arkode_mem, ARKLsMassFn mass)

   Specifies the mass matrix approximation routine to be used for the
   matrix-based solver with the ARKLS interface.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mass: name of user-supplied mass matrix approximation function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_MASSMEM_NULL: the mass matrix solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an argument had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This routine must be called after the ARKLS mass matrix
      solver interface has been initialized through a call to
      :c:func:`ARKodeSetMassLinearSolver`.

      Since there is no default difference quotient function for mass
      matrices, *mass* must be non-``NULL``.

      The function type :c:func:`ARKLsMassFn` is described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff)

   Enables or disables scaling the linear system solution to account for a
   change in :math:`\gamma` in the linear system. For more details see
   :numref:`SUNLinSol.Lagged_matrix`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param onoff: flag to enable (``SUNTRUE``) or disable (``SUNFALSE``)
                 scaling.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_ILL_INPUT: the attached linear solver is not matrix-based.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Linear solution scaling is enabled by default when a matrix-based
      linear solver is attached.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.linear_solution_scaling".

   .. versionadded:: 6.1.0


.. _ARKODE.Usage.ARKLsInputs.MatrixFree:

Optional inputs for matrix-free ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==================================================  =================================  ==================
Optional input                                      Function name                      Default
==================================================  =================================  ==================
:math:`Jv` functions (*jtimes* and *jtsetup*)       :c:func:`ARKodeSetJacTimes`        DQ,  none
:math:`Jv` DQ rhs function (*jtimesRhsFn*)          :c:func:`ARKodeSetJacTimesRhsFn`   fi
:math:`Mv` functions (*mtimes* and *mtsetup*)       :c:func:`ARKodeSetMassTimes`       none, none
==================================================  =================================  ==================


As described in :numref:`ARKODE.Mathematics.Linear`, when solving
the Newton linear systems with matrix-free methods, the ARKLS
interface requires a *jtimes* function to compute an approximation to
the product between the Jacobian matrix
:math:`J(t,y)` and a vector :math:`v`. The user can supply a custom
Jacobian-times-vector approximation function, or use the default
internal difference quotient function that comes with the ARKLS
interface.

A user-defined Jacobian-vector function must be of type
:c:type:`ARKLsJacTimesVecFn` and can be specified through a call
to :c:func:`ARKodeSetJacTimes` (see :numref:`ARKODE.Usage.UserSupplied`
for specification details).  As with the
user-supplied preconditioner functions, the evaluation and
processing of any Jacobian-related data needed by the user's
Jacobian-times-vector function is done in the optional user-supplied
function of type :c:type:`ARKLsJacTimesSetupFn` (see
:numref:`ARKODE.Usage.UserSupplied` for specification details).  As with
the preconditioner functions, a pointer to the user-defined
data structure, *user_data*, specified through
:c:func:`ARKodeSetUserData` (or a ``NULL`` pointer otherwise) is
passed to the Jacobian-times-vector setup and product functions each
time they are called.


.. c:function:: int ARKodeSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup, ARKLsJacTimesVecFn jtimes)

   Specifies the Jacobian-times-vector setup and product functions.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param jtsetup: user-defined Jacobian-vector setup function.
                   Pass ``NULL`` if no setup is necessary.
   :param jtimes: user-defined Jacobian-vector product function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARKLS_SUNLS_FAIL: an error occurred when setting up
                             the Jacobian-vector product in the ``SUNLinearSolver``
                             object used by the ARKLS interface.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default is to use an internal finite difference
      quotient for *jtimes* and to leave out *jtsetup*.  If ``NULL`` is
      passed to *jtimes*, these defaults are used.  A user may
      specify non-``NULL`` *jtimes* and ``NULL`` *jtsetup* inputs.

      This function must be called *after* the ARKLS system solver
      interface has been initialized through a call to
      :c:func:`ARKodeSetLinearSolver`.

      The function types :c:type:`ARKLsJacTimesSetupFn` and
      :c:type:`ARKLsJacTimesVecFn` are described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0


When using the internal difference quotient the user may optionally supply
an alternative implicit right-hand side function for use in the Jacobian-vector
product approximation by calling :c:func:`ARKodeSetJacTimesRhsFn`. The
alternative implicit right-hand side function should compute a suitable (and
differentiable) approximation to the :math:`f^I` function provided to
``*StepCreate``. For example, as done in :cite:p:`dorr2010numerical`,
the alternative function may use lagged values when evaluating a nonlinearity
in :math:`f^I` to avoid differencing a potentially non-differentiable factor.
We note that in many instances this same :math:`f^I` routine would also have
been desirable for the nonlinear solver, in which case the user should specify
this through calls to *both* :c:func:`ARKodeSetJacTimesRhsFn` and
:c:func:`ARKodeSetNlsRhsFn`.


.. c:function:: int ARKodeSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn)

   Specifies an alternative implicit right-hand side function for use in the
   internal Jacobian-vector product difference quotient approximation.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param jtimesRhsFn: the name of the C function (of type
                       :c:func:`ARKRhsFn`) defining the alternative right-hand side function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default is to use the implicit right-hand side function
      provided to ``*StepCreate`` in the internal difference quotient. If
      the input implicit right-hand side function is ``NULL``, the default is used.

      This function must be called *after* the ARKLS system solver interface has
      been initialized through a call to :c:func:`ARKodeSetLinearSolver`.

   .. versionadded:: 6.1.0


Similarly, if a problem involves a non-identity mass matrix,
:math:`M\ne I`, then matrix-free solvers require a *mtimes* function
to compute an approximation to the product between the mass matrix
:math:`M(t)` and a vector :math:`v`.  This function must be
user-supplied since there is no default value, it must be
of type :c:func:`ARKLsMassTimesVecFn`, and can be specified
through a call to the  :c:func:`ARKodeSetMassTimes` routine.
Similarly to the user-supplied preconditioner functions, any evaluation
and processing of any mass matrix-related data needed by the user's
mass-matrix-times-vector function may be done in an optional user-supplied
function of type :c:type:`ARKLsMassTimesSetupFn` (see
:numref:`ARKODE.Usage.UserSupplied` for specification details).



.. c:function:: int ARKodeSetMassTimes(void* arkode_mem, ARKLsMassTimesSetupFn mtsetup, ARKLsMassTimesVecFn mtimes, void* mtimes_data)

   Specifies the mass matrix-times-vector setup and product functions.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mtsetup: user-defined mass matrix-vector setup function.
                   Pass ``NULL`` if no setup is necessary.
   :param mtimes: user-defined mass matrix-vector product function.
   :param mtimes_data: a pointer to user data, that will be supplied
                       to both the *mtsetup* and *mtimes* functions.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_MASSMEM_NULL: the mass matrix solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARKLS_SUNLS_FAIL: an error occurred when setting up
                             the mass-matrix-vector product in the ``SUNLinearSolver``
                             object used by the ARKLS interface.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      There is no default finite difference quotient for
      *mtimes*, so if using the ARKLS mass matrix solver interface with
      NULL-valued SUNMATRIX input :math:`M`, and this routine is called
      with NULL-valued *mtimes*, an error will occur.  A user may
      specify ``NULL`` for *mtsetup*.

      This function must be called *after* the ARKLS mass
      matrix solver interface has been initialized through a call to
      :c:func:`ARKodeSetMassLinearSolver`.

      The function types :c:type:`ARKLsMassTimesSetupFn` and
      :c:type:`ARKLsMassTimesVecFn` are described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.ARKLsInputs.Iterative:

Optional inputs for iterative ``SUNLinearSolver`` modules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

====================================================  ======================================  ==================
Optional input                                        Function name                           Default
====================================================  ======================================  ==================
Newton preconditioning functions                      :c:func:`ARKodeSetPreconditioner`       ``NULL``, ``NULL``
Mass matrix preconditioning functions                 :c:func:`ARKodeSetMassPreconditioner`   ``NULL``, ``NULL``
Newton linear and nonlinear tolerance ratio           :c:func:`ARKodeSetEpsLin`               0.05
Mass matrix linear and nonlinear tolerance ratio      :c:func:`ARKodeSetMassEpsLin`           0.05
Newton linear solve tolerance conversion factor       :c:func:`ARKodeSetLSNormFactor`         vector length
Mass matrix linear solve tolerance conversion factor  :c:func:`ARKodeSetMassLSNormFactor`     vector length
====================================================  ======================================  ==================


As described in :numref:`ARKODE.Mathematics.Linear`, when using
an iterative linear solver the user may supply a preconditioning
operator to aid in solution of the system.  This operator consists of
two user-supplied functions, *psetup* and *psolve*, that are supplied
to ARKODE using either the function
:c:func:`ARKodeSetPreconditioner` (for preconditioning the
Newton system), or the function
:c:func:`ARKodeSetMassPreconditioner` (for preconditioning the
mass matrix system).  The *psetup* function supplied to these routines
should handle evaluation and preprocessing of any Jacobian or
mass-matrix data needed by the user's preconditioner solve function,
*psolve*.  The user data pointer received through
:c:func:`ARKodeSetUserData` (or a pointer to ``NULL`` if user data
was not specified) is passed to the *psetup* and *psolve* functions.
This allows the user to create an arbitrary
structure with relevant problem data and access it during the
execution of the user-supplied preconditioner functions without using
global data in the program.  If preconditioning is supplied for both
the Newton and mass matrix linear systems, it is expected that the
user will supply different *psetup* and *psolve* function for each.

Also, as described in :numref:`ARKODE.Mathematics.Error.Linear`, the
ARKLS interface requires that iterative linear solvers stop when
the norm of the preconditioned residual satisfies

.. math::
   \|r\| \le \frac{\epsilon_L \epsilon}{10}

where the default :math:`\epsilon_L = 0.05` may be modified by
the user through the :c:func:`ARKodeSetEpsLin` function.


.. c:function:: int ARKodeSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup, ARKLsPrecSolveFn psolve)

   Specifies the user-supplied preconditioner setup and solve functions.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param psetup: user defined preconditioner setup function.  Pass
                  ``NULL`` if no setup is needed.
   :param psolve: user-defined preconditioner solve function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARKLS_SUNLS_FAIL: an error occurred when setting up preconditioning
                             in the ``SUNLinearSolver`` object used
                             by the ARKLS interface.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The default is ``NULL`` for both arguments (i.e., no
      preconditioning).

      This function must be called *after* the ARKLS system solver
      interface has been initialized through a call to
      :c:func:`ARKodeSetLinearSolver`.

      Both of the function types :c:func:`ARKLsPrecSetupFn` and
      :c:func:`ARKLsPrecSolveFn` are described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMassPreconditioner(void* arkode_mem, ARKLsMassPrecSetupFn psetup, ARKLsMassPrecSolveFn psolve)

   Specifies the mass matrix preconditioner setup and solve functions.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param psetup: user defined preconditioner setup function.  Pass
                  ``NULL`` if no setup is to be done.
   :param psolve: user-defined preconditioner solve function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARKLS_SUNLS_FAIL: an error occurred when setting up preconditioning
                             in the ``SUNLinearSolver`` object used
                             by the ARKLS interface.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This function must be called *after* the ARKLS mass
      matrix solver interface has been initialized through a call to
      :c:func:`ARKodeSetMassLinearSolver`.

      The default is ``NULL`` for both arguments (i.e. no
      preconditioning).

      Both of the function types :c:func:`ARKLsMassPrecSetupFn` and
      :c:func:`ARKLsMassPrecSolveFn` are described in
      :numref:`ARKODE.Usage.UserSupplied`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetEpsLin(void* arkode_mem, sunrealtype eplifac)

   Specifies the factor :math:`\epsilon_L` by which the tolerance on
   the nonlinear iteration is multiplied to get a tolerance on the
   linear iteration.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param eplifac: linear convergence safety factor.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      Passing a value *eplifac* :math:`\le 0` indicates to use the
      default value of 0.05.

      This function must be called *after* the ARKLS system solver
      interface has been initialized through a call to
      :c:func:`ARKodeSetLinearSolver`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.eps_lin".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMassEpsLin(void* arkode_mem, sunrealtype eplifac)

   Specifies the factor by which the tolerance on the nonlinear
   iteration is multiplied to get a tolerance on the mass matrix
   linear iteration.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param eplifac: linear convergence safety factor.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_MASSMEM_NULL: the mass matrix solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This function must be called *after* the ARKLS mass
      matrix solver interface has been initialized through a call to
      :c:func:`ARKodeSetMassLinearSolver`.

      Passing a value *eplifac* :math:`\le 0` indicates to use the default value
      of 0.05.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.mass_eps_lin".

   .. versionadded:: 6.1.0


Since iterative linear solver libraries typically consider linear residual
tolerances using the :math:`L_2` norm, whereas ARKODE focuses on errors
measured in the WRMS norm :eq:`ARKODE_WRMS_NORM`, the ARKLS interface internally
converts between these quantities when interfacing with linear solvers,

.. math::
   \text{tol}_{L2} = \textit{nrmfac}\ \ \text{tol}_{WRMS}.
   :label: ARKODE_NRMFAC

Prior to the introduction of :c:func:`N_VGetLength` in SUNDIALS v5.0.0 the
value of :math:`nrmfac` was computed using the vector dot product.  Now, the
functions :c:func:`ARKodeSetLSNormFactor` and :c:func:`ARKodeSetMassLSNormFactor`
allow for additional user control over these conversion factors.


.. c:function:: int ARKodeSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac)

   Specifies the factor to use when converting from the integrator tolerance
   (WRMS norm) to the linear solver tolerance (L2 norm) for Newton linear system
   solves.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nrmfac: the norm conversion factor. If *nrmfac* is:

                  :math:`> 0` then the provided value is used.

                  :math:`= 0` then the conversion factor is computed using the vector
                  length i.e., ``nrmfac = sqrt(N_VGetLength(y))`` (*default*).

                  :math:`< 0` then the conversion factor is computed using the vector dot
                  product i.e., ``nrmfac = sqrt(N_VDotProd(v,v))`` where all the entries
                  of ``v`` are one.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported by the
                                    current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This function must be called *after* the ARKLS system solver interface has
      been initialized through a call to :c:func:`ARKodeSetLinearSolver`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.ls_norm_factor".

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetMassLSNormFactor(void* arkode_mem, sunrealtype nrmfac)

   Specifies the factor to use when converting from the integrator tolerance
   (WRMS norm) to the linear solver tolerance (L2 norm) for mass matrix linear
   system solves.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nrmfac: the norm conversion factor. If *nrmfac* is:

                  :math:`> 0` then the provided value is used.

                  :math:`= 0` then the conversion factor is computed using the vector
                  length i.e., ``nrmfac = sqrt(N_VGetLength(y))`` (*default*).

                  :math:`< 0` then the conversion factor is computed using the vector dot
                  product i.e., ``nrmfac = sqrt(N_VDotProd(v,v))`` where all the entries
                  of ``v`` are one.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This function must be called *after* the ARKLS mass matrix solver interface
      has been initialized through a call to :c:func:`ARKodeSetMassLinearSolver`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.mass_ls_norm_factor".

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.ARKodeRootfindingInputTable:


Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to
control the rootfinding algorithm, the mathematics of which are
described in :numref:`ARKODE.Mathematics.Rootfinding`.


.. cssclass:: table-bordered

======================================  =====================================  ==================
Optional input                          Function name                          Default
======================================  =====================================  ==================
Direction of zero-crossings to monitor  :c:func:`ARKodeSetRootDirection`       both
Disable inactive root warnings          :c:func:`ARKodeSetNoInactiveRootWarn`  enabled
======================================  =====================================  ==================



.. c:function:: int ARKodeSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rootdir: state array of length *nrtfn*, the number of root
                   functions :math:`g_i` (the value of *nrtfn* was supplied in
                   the call to :c:func:`ARKodeRootInit`).  If ``rootdir[i] ==
                   0`` then crossing in either direction for :math:`g_i` should be
                   reported.  A value of +1 or -1 indicates that the solver
                   should report only zero-crossings where :math:`g_i` is
                   increasing or decreasing, respectively.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      The default behavior is to monitor for both zero-crossing directions.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.

   :param arkode_mem: pointer to the ARKODE memory block.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      ARKODE will not report the initial conditions as a
      possible zero-crossing (assuming that one or more components
      :math:`g_i` are zero at the initial time).  However, if it appears
      that some :math:`g_i` is identically zero at the initial time
      (i.e., :math:`g_i` is zero at the initial time *and* after the
      first step), ARKODE will issue a warning which can be disabled with
      this optional input function.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.no_inactive_root_warn".

   .. versionadded:: 6.1.0




.. _ARKODE.Usage.ARKodeProcessingInputTable:

Step and stage processing optional inputs (ADVANCED)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ARKODE provides multiple options for user-supplied callback routines
that can be called at various times within the time-stepping process.
Each of these callback functions has the same structure, wherein the
callback function will be provided with the current time, current
solution, and the *user_data* structure that was provided to
:c:func:`ARKodeSetUserData`.  Specifically, users may provide
separate callback functions for the following events within a time
step:

* just prior to starting a time step attempt
  (:c:func:`ARKodeSetPreprocessStepFn`),

* at the end of a successful time step, but before overwriting the
  "saved" state (:c:func:`ARKodeSetPostprocessStepFn`),

* at the end of a failed time step (e.g., but before the step is
  reattempted (:c:func:`ARKodeSetPostprocessStepFailFn`),

* just prior to evaluating user-provided right-hand side (RHS)
  functions (:c:func:`ARKodeSetPreprocessRHSFn`)

* immediately after each stage is completed within a time step
  (:c:func:`ARKodeSetPostprocessStageFn`)

For users who wish to perform different actions at individual internal
stages within an ARKODE method, they may obtain the current stage index by
calling :c:func:`ARKodeGetStageIndex`
in their stage-level callback routines provided to
:c:func:`ARKodeSetPreprocessRHSFn` and :c:func:`ARKodeSetPostprocessStageFn`.

The specific ordering of these functions within a given step depends on
whether each stage is explicit (as in ERKStep) or implicit (as in ARKStep or
MRIStep).  Denoting the last "saved" time step as :math:`(t_n,y_n)`, the
time-evolving state within a step as :math:`(t_{cur},y_{cur})`, the functions
provided to the five above functions as ``PreprocessStep``, ``PostprocessStep``,
``PostprocessFailedStep``, ``PreprocessRHS``, and ``PostprocessStage``, and
denoting the IVP right hand side function as ``RHS``, then
the flow of a 3-stage explicit method would proceed as:

1. Update :math:`(t_{cur},y_{cur}) \gets (t_n,y_n)`
2. Call ``PreprocessStep`` with :math:`(t_{cur},y_{cur})`
3. Stage 0
   a. Call ``PreprocessRHS`` with :math:`(t_{cur},y_{cur})`
   b. Evaluate ``RHS`` at :math:`(t_{cur},y_{cur})`
   c. Update :math:`(t_{cur},y_{cur})` with the next stage solution
   d. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
4. Stage 1
   a. Call ``PreprocessRHS`` with :math:`(t_{cur},y_{cur})`
   b. Evaluate ``RHS`` at :math:`(t_{cur},y_{cur})`
   c. Update :math:`(t_{cur},y_{cur})` with the next stage solution
   d. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
5. Stage 2
   a. Call ``PreprocessRHS`` with :math:`(t_{cur},y_{cur})`
   b. Evaluate ``RHS`` at :math:`(t_{cur},y_{cur})`
   c. Update :math:`(t_{cur},y_{cur})` with the new time step solution
   d. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
6. Check the local error.
   a. If the step is successful then call ``PostProcessStep`` with  :math:`(t_{cur},y_{cur})`, determine the next internal step size :math:`h_n`, and update :math:`(t_n,y_n) \gets (t_{cur},y_{cur})`
   b. Else call ``PostProcessFailedStep`` with  :math:`(t_{cur},y_{cur})`, determine the next internal step size :math:`h_n`, and return to step 3.


Alternately, the flow of a 3-stage method that must perform a solve of some sort for each stage (i.e., a DIRK or ARK method in ARKStep, or a multirate method with MRIstep) would proceed as follows.  Here, we show the implicit-explicit approach since that also shows the relationship between both the implicit right-hand side function ``RHS_i`` and the explicit right-hand side function ``RHS_e``:

1. Update :math:`(t_{cur},y_{cur}) \gets (t_n,y_n)`
2. Call ``PreprocessStep`` with :math:`(t_{cur},y_{cur})`
3. Stage 0
   a. Solve implicit system, calling ``PreprocessRHS`` and then ``RHS_i`` with :math:`(t_{cur},y_{cur})` at each solver iteration; at the end of this iteration :math:`(t_{cur},y_{cur})` holds the updated stage solution
   b. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
   c. Call ``PreprocessRHS`` with :math:`(t_{cur},y_{cur})`
   d. Evaluate ``RHS_i`` and then ``RHS_e`` at :math:`(t_{cur},y_{cur})`
4. Stage 1
   a. Solve implicit system, calling ``PreprocessRHS`` and then ``RHS_i`` with :math:`(t_{cur},y_{cur})` at each solver iteration; at the end of this iteration :math:`(t_{cur},y_{cur})` holds the updated stage solution
   b. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
   c. Call ``PreprocessRHS`` with :math:`(t_{cur},y_{cur})`
   d. Evaluate ``RHS_i`` and then ``RHS_e`` at :math:`(t_{cur},y_{cur})`
5. Stage 2
   a. Solve implicit system, calling ``PreprocessRHS`` and then ``RHS_i`` with :math:`(t_{cur},y_{cur})` at each solver iteration; at the end of this iteration :math:`(t_{cur},y_{cur})` holds the updated stage solution
   b. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
   c. Call ``PreprocessRHS`` with :math:`(t_{cur},y_{cur})`
   d. Evaluate ``RHS_i`` and then ``RHS_e`` at :math:`(t_{cur},y_{cur})`
6. Update :math:`(t_{cur},y_{cur})` with the new time step solution
7. Call ``PostprocessStage`` with :math:`(t_{cur},y_{cur})`
8. Check the local error.
   a. If the step is successful then call ``PostProcessStep`` with  :math:`(t_{cur},y_{cur})`, determine the next internal step size :math:`h_n`, and update :math:`(t_n,y_n) \gets (t_{cur},y_{cur})`
   b. Else call ``PostProcessFailedStep`` with  :math:`(t_{cur},y_{cur})`, determine the next internal step size :math:`h_n`, and return to step 3.

We consider these as "advanced" because of their danger, although the
callback functions are provided with the internally-evolving state,
users should **not** adjust entries of this state vector, since doing
so will destroy all theoretical guarantees of solution accuracy and
numerical stability.


.. cssclass:: table-bordered

=================================================  ==========================================  =======================
Optional input                                     Function name                               Default
=================================================  ==========================================  =======================
Set time step preprocessing function               :c:func:`ARKodeSetPreprocessStepFn`         ``NULL``
Set time step postprocessing function              :c:func:`ARKodeSetPostprocessStepFn`        ``NULL``
Set failed time step postrocessing function        :c:func:`ARKodeSetPostprocessStepFailFn`    ``NULL``
Set right-hand side preprocessing function         :c:func:`ARKodeSetPreprocessRHSFn`          ``NULL``
Set stage postprocessing function.                 :c:func:`ARKodeSetPostprocessStageFn`       ``NULL``
=================================================  ==========================================  =======================



.. c:function:: int ARKodeSetPreprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)

   [ADVANCED] Provides a function to be called before each time step.  A ``NULL``
   input function disables step preprocessing.  If a user-supplied failed-step postprocessing
   function is supplied by calling :c:func:`ARKodeSetPostprocessStepFailFn`, then
   preprocessing will not be called on the step attempt that immediately follows a failed step.

   This should **not** adjust the state vector itself.  It is designed to allow users to set
   up auxiliary data structures that they will use within the time step (e.g., in their
   right-hand side function(s)).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ProcessStep: the user-supplied function to call.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. danger::

      If the supplied function modifies any of the active state data, then all
      theoretical guarantees of solution accuracy and stability are lost.

   .. versionadded:: x.y.z


.. c:function:: int ARKodeSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)

   [ADVANCED] Provides a function to be called following each successful time step.  A ``NULL``
   input function disables step postprocessing.

   This should **not** adjust the state vector itself.  It is designed to allow users to compute
   relevant diagnostic information after each step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ProcessStep: the user-supplied function to call.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. danger::

      If the supplied function modifies any of the active state data, then all
      theoretical guarantees of solution accuracy and stability are lost.


.. c:function:: int ARKodeSetPostprocessStepFailFn(void* arkode_mem, ARKPostProcessFn ProcessStep)

   [ADVANCED] Provides a function to be called following each failed time step.  A ``NULL``
   input function disables failed step postprocessing.

   This should **not** adjust the state vector itself.  It is designed to allow users to reset
   any relevant diagnostic information they may have accumulated within a rejected time step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ProcessStep: the user-supplied function to call.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. danger::

      If the supplied function modifies any of the active state data, then all
      theoretical guarantees of solution accuracy and stability are lost.

   .. versionadded:: x.y.z


.. c:function:: int ARKodeSetPreprocessRHSFn(void* arkode_mem, ARKPostProcessFn ProcessStep)

   [ADVANCED] Provides a function to be called prior to evaluating user-provided right-hand
   side (RHS) functions.  For partitioned methods (e.g., ARKStep or MRIStep), that will call
   multiple RHS functions with identical inputs, this is called only once prior to the first
   RHS evaluation.  A ``NULL`` input function disables RHS preprocessing.

   This should **not** adjust the state vector itself.  It is designed to allow users to set up
   auxiliary data structures that will be used within the RHS evaluations (e.g., MPI communication
   to fill and send exchange buffers).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ProcessStep: the user-supplied function to call.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. danger::

      If the supplied function modifies any of the active state data, then all
      theoretical guarantees of solution accuracy and stability are lost.

   .. versionadded:: x.y.z


.. c:function:: int ARKodeSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStep)

   [ADVANCED] Provides a function to be called immediately after each stage is completed within
   ARKODE's multi-stage methods.  A ``NULL`` input function disables stage postprocessing.

   This should **not** adjust the state vector itself.  It is designed to allow users to compute
   relevant diagnostic information within each step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ProcessStep: the user-supplied function to call.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. danger::

      If the supplied function modifies any of the active state data, then all
      theoretical guarantees of solution accuracy and stability are lost.





.. _ARKODE.Usage.InterpolatedOutput:

Interpolated output function
--------------------------------

An optional function :c:func:`ARKodeGetDky` is available to obtain
additional values of solution-related quantities.  This function
should only be called after a successful return from
:c:func:`ARKodeEvolve`, as it provides interpolated values either of
:math:`y` or of its derivatives (up to the 5th derivative)
interpolated to any value of :math:`t` in the last internal step taken
by :c:func:`ARKodeEvolve`.  Internally, this "dense output" or
"continuous extension" algorithm is identical to the algorithm used for
the maximum order implicit predictors, described in
:numref:`ARKODE.Mathematics.Predictors.Max`, except that derivatives of the
polynomial model may be evaluated upon request.



.. c:function:: int ARKodeGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function
   :math:`y` at the time *t*,
   i.e. :math:`y^{(k)}(t)`, for values of the
   independent variable satisfying :math:`t_n-h_n \le t \le t_n`, with
   :math:`t_n` as current internal time reached, and :math:`h_n` is
   the last internal step size successfully used by the solver.  This
   routine uses an interpolating polynomial of degree *min(degree, 5)*,
   where *degree* is the argument provided to
   :c:func:`ARKodeSetInterpolantDegree`.  The user may request *k* in the
   range {0,..., *min(degree, kmax)*} where *kmax* depends on the choice of
   interpolation module. For Hermite interpolants *kmax = 5* and for Lagrange
   interpolants *kmax = 3*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param t: the value of the independent variable at which the
             derivative is to be evaluated.
   :param k: the derivative order requested.
   :param dky: output vector (must be allocated by the user).

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_BAD_K: *k* is not in the range {0,..., *min(degree, kmax)*}.
   :retval ARK_BAD_T: *t* is not in the interval :math:`[t_n-h_n, t_n]`.
   :retval ARK_BAD_DKY: the *dky* vector was ``NULL``.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      It is only legal to call this function after a successful
      return from :c:func:`ARKodeEvolve`.

      A user may access the values :math:`t_n` and :math:`h_n` via the
      functions :c:func:`ARKodeGetCurrentTime` and
      :c:func:`ARKodeGetLastStep`, respectively.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.OptionalOutputs:

Optional output functions
------------------------------

ARKODE provides an extensive set of functions that can be used to
obtain solver performance information.  We organize these into groups:

#. General ARKODE output routines are in
   :numref:`ARKODE.Usage.ARKodeMainOutputs`,
#. ARKODE implicit solver output routines are in
   :numref:`ARKODE.Usage.ARKodeImplicitSolverOutputs`,
#. Output routines regarding root-finding results are in
   :numref:`ARKODE.Usage.ARKodeRootOutputs`,
#. Linear solver output routines are in
   :numref:`ARKODE.Usage.ARKLsOutputs` and
#. General usability routines (e.g. to print the current ARKODE
   parameters, or output the current Butcher table(s)) are in
   :numref:`ARKODE.Usage.ARKodeExtraOutputs`.

Following each table, we elaborate on each function.

Some of the optional outputs, especially the various counters, can be
very useful in determining the efficiency of various methods inside
ARKODE.  For example:

* The counters *nsteps*, *nfe_evals* and *nfi_evals*
  provide a rough measure of the overall cost of a given run, and can
  be compared between runs with different solver options to suggest
  which set of options is the most efficient.

* The ratio *nniters/nsteps* measures the performance of the
  nonlinear iteration in solving the nonlinear systems at each stage,
  providing a measure of the degree of nonlinearity in the problem.
  Typical values of this for a Newton solver on a general problem
  range from 1.1 to 1.8.

* When using a Newton nonlinear solver, the ratio *njevals/nniters*
  (when using a direct linear solver), and the ratio
  *nliters/nniters* (when using an iterative linear solver) can
  indicate the quality of the approximate Jacobian or preconditioner being
  used.  For example, if this ratio is larger for a user-supplied
  Jacobian or Jacobian-vector product routine than for the
  difference-quotient routine, it can indicate that the user-supplied
  Jacobian is inaccurate.

* The ratio *expsteps/accsteps* can measure the quality of the ImEx
  splitting used, since a higher-quality splitting will be dominated
  by accuracy-limited steps, and hence a lower ratio.

* The ratio *nsteps/step_attempts* can measure the quality of the
  time step adaptivity algorithm, since a poor algorithm will result
  in more failed steps, and hence a lower ratio.

It is therefore recommended that users retrieve and output these
statistics following each run, and take some time to investigate
alternate solver options that will be more optimal for their
particular problem of interest.



.. _ARKODE.Usage.ARKodeMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=====================================================  ============================================
Optional output                                        Function name
=====================================================  ============================================
Size of ARKODE real and integer workspaces             :c:func:`ARKodeGetWorkSpace`
Cumulative number of internal steps                    :c:func:`ARKodeGetNumSteps`
Actual initial time step size used                     :c:func:`ARKodeGetActualInitStep`
Step size used for the last successful step            :c:func:`ARKodeGetLastStep`
Step size to be attempted on the next step             :c:func:`ARKodeGetCurrentStep`
Integration direction, e.g., forward or backward       :c:func:`ARKodeGetStepDirection`
Last saved time reached by the solver                  :c:func:`ARKodeGetLastTime`
Last saved solution reached by the solver              :c:func:`ARKodeGetLastState`
Current internal time reached by the solver            :c:func:`ARKodeGetCurrentTime`
Current internal solution reached by the solver        :c:func:`ARKodeGetCurrentState`
Current :math:`\gamma` value used by the solver        :c:func:`ARKodeGetCurrentGamma`
Suggested factor for tolerance scaling                 :c:func:`ARKodeGetTolScaleFactor`
Error weight vector for state variables                :c:func:`ARKodeGetErrWeights`
Residual weight vector                                 :c:func:`ARKodeGetResWeights`
Single accessor to many statistics at once             :c:func:`ARKodeGetStepStats`
Print all statistics                                   :c:func:`ARKodePrintAllStats`
Name of constant associated with a return flag         :c:func:`ARKodeGetReturnFlagName`
No. of explicit stability-limited steps                :c:func:`ARKodeGetNumExpSteps`
No. of accuracy-limited steps                          :c:func:`ARKodeGetNumAccSteps`
No. of attempted steps                                 :c:func:`ARKodeGetNumStepAttempts`
No. of RHS evaluations                                 :c:func:`ARKodeGetNumRhsEvals`
No. of local error test failures that have occurred    :c:func:`ARKodeGetNumErrTestFails`
No. of failed steps due to a nonlinear solver failure  :c:func:`ARKodeGetNumStepSolveFails`
Estimated local truncation error vector                :c:func:`ARKodeGetEstLocalErrors`
Number of constraint test failures                     :c:func:`ARKodeGetNumConstrFails`
Retrieve a pointer for user data                       :c:func:`ARKodeGetUserData`
Retrieve the accumulated temporal error estimate       :c:func:`ARKodeGetAccumulatedError`
Current stage index, and total number of stages        :c:func:`ARKodeGetStageIndex`
=====================================================  ============================================




.. c:function:: int ARKodeGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)

   Returns the ARKODE real and integer workspace sizes.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lenrw: the number of ``sunrealtype`` values in the ARKODE workspace.
   :param leniw: the number of integer values in the ARKODE workspace.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.


.. c:function:: int ARKodeGetNumSteps(void* arkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nsteps: number of steps taken in the solver.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetActualInitStep(void* arkode_mem, sunrealtype* hinused)

   Returns the value of the integration step size used on the first step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hinused: actual value of initial step size.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      Even if the value of the initial integration step was
      specified by the user through a call to
      :c:func:`ARKodeSetInitStep`, this value may have been changed by
      ARKODE to ensure that the step size fell within the prescribed
      bounds :math:`(h_{min} \le h_0 \le h_{max})`, or to satisfy the
      local error test condition, or to ensure convergence of the
      nonlinear solver.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetLastStep(void* arkode_mem, sunrealtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hlast: step size taken on the last internal step.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetCurrentStep(void* arkode_mem, sunrealtype* hcur)

   Returns the integration step size to be attempted on the next internal step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param hcur: step size to be attempted on the next internal step.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetStepDirection(void *arkode_mem, sunrealtype *stepdir)

   Returns the direction of integration that will be used on the next internal
   step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param stepdir: a positive number if integrating forward, a negative number
                   if integrating backward, or zero if the direction has not
                   been set.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.2.0


.. c:function:: int ARKodeGetLastTime(void* arkode_mem, sunrealtype* tn)

   Returns the last saved time reached by the solver.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param tn: last saved time reached.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: x.y.z


.. c:function:: int ARKodeGetLastState(void *arkode_mem, N_Vector *yn)

   Returns the last saved solution reached by the solver.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param yn: last saved solution.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      Users should exercise extreme caution when using this function,
      as altering values of *yn* may lead to undesirable behavior, depending
      on the particular use case and on when this routine is called.

   .. versionadded:: x.y.z


.. c:function:: int ARKodeGetCurrentTime(void* arkode_mem, sunrealtype* tcur)

   Returns the current internal time reached by the solver.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param tcur: current internal time reached.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetCurrentState(void *arkode_mem, N_Vector *ycur)

   Returns the current internal solution reached by the solver.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ycur: current internal solution.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      Users should exercise extreme caution when using this function,
      as altering values of *ycur* may lead to undesirable behavior, depending
      on the particular use case and on when this routine is called.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetCurrentGamma(void *arkode_mem, sunrealtype *gamma)

   Returns the current internal value of :math:`\gamma` used in the implicit
   solver Newton matrix (see equation :eq:`ARKODE_NewtonMatrix`).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param gamma: current step size scaling factor in the Newton system.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param tolsfac: suggested scaling factor for user-supplied tolerances.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetErrWeights(void* arkode_mem, N_Vector eweight)

   Returns the current error weight vector.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param eweight: solution error weights at the current time.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      The user must allocate space for *eweight*, that will be
      filled in by this function.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetResWeights(void* arkode_mem, N_Vector rweight)

   Returns the current residual weight vector.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rweight: residual error weights at the current time.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      The user must allocate space for *rweight*, that will be
      filled in by this function.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused, sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)

   Returns many of the most useful optional outputs in a single call.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nsteps: number of steps taken in the solver.
   :param hinused: actual value of initial step size.
   :param hlast: step size taken on the last internal step.
   :param hcur: step size to be attempted on the next internal step.
   :param tcur: current internal time reached.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodePrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)

   Outputs all of the integrator, nonlinear solver, linear solver, and other
   statistics.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param outfile: pointer to output file.
   :param fmt: the output format:

               * :c:enumerator:`SUN_OUTPUTFORMAT_TABLE` -- prints a table of values

               * :c:enumerator:`SUN_OUTPUTFORMAT_CSV` -- prints a comma-separated list
                 of key and value pairs e.g., ``key1,value1,key2,value2,...``

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: an invalid formatting option was provided.

   .. note::

      The Python module ``tools/suntools`` provides utilities to read and output
      the data from a SUNDIALS CSV output file using the key and value pair
      format.

   .. versionadded:: 6.1.0


.. c:function:: char* ARKodeGetReturnFlagName(long int flag)

   Returns the name of the ARKODE constant corresponding to *flag*.
   See :ref:`ARKODE.Constants`.

   :param flag: a return flag from an ARKODE function.

   :return: The return value is a string containing the name of
            the corresponding constant.

   .. warning::

      The user is responsible for freeing the returned string.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumExpSteps(void* arkode_mem, long int* expsteps)

   Returns the cumulative number of stability-limited steps
   taken by the solver (so far). If the combination of the maximum number of stages
   and the current time step size in the LSRKStep module will not allow for a stable
   step, the counter also accounts for such returns.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param expsteps: number of stability-limited steps taken in the solver.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumAccSteps(void* arkode_mem, long int* accsteps)

   Returns the cumulative number of accuracy-limited steps
   taken by the solver (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param accsteps: number of accuracy-limited steps taken in the solver.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param step_attempts: number of steps attempted by solver.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumRhsEvals(void* arkode_mem, int partition_index, long int* num_rhs_evals)

   Returns the number of calls to the user's right-hand side function (so far).
   For implicit methods or methods with an implicit partition, the count does
   not include calls made by a linear solver or preconditioner.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param num_partition: the right-hand side partition index:

                   * For ERKStep, ``0`` corresponds to :math:`f(t,y)`

                   * For ARKStep, ``0`` corresponds to :math:`f^E(t,y)` and
                     ``1`` to :math:`f^I(t,y)`

                   * For MRIStep, ``0`` corresponds to :math:`f^E(t,y)` and
                     ``1`` to :math:`f^I(t,y)`

                   * For SPRKStep, ``0`` corresponds to :math:`f_1(t,p)` and
                     ``1`` to :math:`f_2(t,q)`

                   A negative index will return the sum of the evaluations for
                   each partition.

   :param num_rhs_evals: the number of right-hand side evaluations.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: if ``arkode_mem`` was ``NULL``.
   :retval ARK_ILL_INPUT: if ``num_partiton`` was invalid for the stepper or
                          ``num_rhs_evals`` was ``NULL``

   .. versionadded:: 6.2.0


.. c:function:: int ARKodeGetNumErrTestFails(void* arkode_mem, long int* netfails)

   Returns the number of local error test failures that
   have occurred (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param netfails: number of error test failures.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumStepSolveFails(void* arkode_mem, long int* ncnf)

   Returns the number of failed steps due to a nonlinear solver failure (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ncnf: number of step failures.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: implicit solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetEstLocalErrors(void* arkode_mem, N_Vector ele)

   Returns the vector of estimated local truncation errors
   for the current step.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ele: vector of estimated local truncation errors.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      The user must allocate space for *ele*, that will be
      filled in by this function.

      The values returned in *ele* are valid only after a successful call
      to :c:func:`ARKodeEvolve` (i.e., it returned a non-negative value).

      The *ele* vector, together with the *eweight* vector from
      :c:func:`ARKodeGetErrWeights`, can be used to determine how the
      various components of the system contributed to the estimated local
      error test.  Specifically, that error test uses the WRMS norm of a
      vector whose components are the products of the components of these
      two vectors.  Thus, for example, if there were recent error test
      failures, the components causing the failures are those with largest
      values for the products, denoted loosely as ``eweight[i]*ele[i]``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumConstrFails(void* arkode_mem, long int* nconstrfails)

   Returns the cumulative number of constraint test failures (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nconstrfails: number of constraint test failures.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: adaptive step sizes are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support temporal adaptivity.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetUserData(void* arkode_mem, void** user_data)

   Returns the user data pointer previously set with
   :c:func:`ARKodeSetUserData`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param user_data: memory reference to a user data pointer.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetAccumulatedError(void* arkode_mem, sunrealtype* accum_error)

   Returns the accumulated temporal error estimate.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param accum_error: pointer to accumulated error estimate.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_WARNING: accumulated error estimation is currently disabled.
   :retval ARK_STEPPER_UNSUPPORTED: temporal error estimation is not supported
                                    by the current time-stepping module.

   .. versionadded:: 6.2.0


.. c:function:: int ARKodeGetStageIndex(void* arkode_mem, int* stage, int *max_stages)

   Returns the index of the current stage (0-based) and the total number of
   stages in the method.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param stage: pointer to storage for the current stage index.
   :param max_stages: pointer to storage for the number of stages in the method.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: x.y.z





.. _ARKODE.Usage.ARKodeImplicitSolverOutputs:

Implicit solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ============================================
Optional output                                      Function name
===================================================  ============================================
Computes state given a correction                    :c:func:`ARKodeComputeState`
Access data to compute the nonlin. sys. function     :c:func:`ARKodeGetNonlinearSystemData`
No. of calls to linear solver setup function         :c:func:`ARKodeGetNumLinSolvSetups`
No. of nonlinear solver iterations                   :c:func:`ARKodeGetNumNonlinSolvIters`
No. of nonlinear solver iterations                   :c:func:`ARKodeGetNumNonlinSolvIters`
No. of nonlinear solver convergence failures         :c:func:`ARKodeGetNumNonlinSolvConvFails`
Single accessor to all nonlinear solver statistics   :c:func:`ARKodeGetNonlinSolvStats`
===================================================  ============================================




.. c:function:: int ARKodeGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)

   Returns the number of calls made to the linear solver's
   setup routine (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nlinsetups: number of linear solver setup calls made.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the nonlinear
      solver object; the counter is reset whenever a new nonlinear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)

   Returns the number of nonlinear solver iterations
   performed (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nniters: number of nonlinear iterations performed.
   :retval ARK_STEPPER_UNSUPPORTED: nonlinear solvers are not supported
                                    by the current time-stepping module.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_NLS_OP_ERR: the SUNNONLINSOL object returned a failure flag.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the nonlinear
      solver object; the counter is reset whenever a new nonlinear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumNonlinSolvConvFails(void* arkode_mem, long int* nncfails)

   Returns the number of nonlinear solver convergence
   failures that have occurred (so far).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nncfails: number of nonlinear convergence failures.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: nonlinear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the nonlinear
      solver object; the counter is reset whenever a new nonlinear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNonlinSolvStats(void* arkode_mem, long int* nniters, long int* nncfails)

   Returns all of the nonlinear solver statistics in a single call.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nniters: number of nonlinear iterations performed.
   :param nncfails: number of nonlinear convergence failures.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARK_NLS_OP_ERR: the SUNNONLINSOL object returned a failure flag.
   :retval ARK_STEPPER_UNSUPPORTED: nonlinear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the nonlinear
      solver object; the counters are reset whenever a new nonlinear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. _ARKODE.Usage.ARKodeRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ==========================================
Optional output                                      Function name
===================================================  ==========================================
Array showing roots found                            :c:func:`ARKodeGetRootInfo`
No. of calls to user root function                   :c:func:`ARKodeGetNumGEvals`
===================================================  ==========================================



.. c:function:: int ARKodeGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to
   have a root.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param rootsfound: array of length *nrtfn* with the indices of the
                      user functions :math:`g_i` found to have a root (the value of
                      *nrtfn* was supplied in the call to
                      :c:func:`ARKodeRootInit`).  For :math:`i = 0 \ldots`
                      *nrtfn*-1, ``rootsfound[i]`` is nonzero if :math:`g_i` has a
                      root, and 0 if not.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      The user must allocate space for *rootsfound* prior to
      calling this function.

      For the components of :math:`g_i` for which a root was found, the
      sign of ``rootsfound[i]`` indicates the direction of
      zero-crossing.  A value of +1 indicates that :math:`g_i` is
      increasing, while a value of -1 indicates a decreasing :math:`g_i`.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ngevals: number of calls made to :math:`g` so far.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.ARKLsOutputs:

Linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A variety of optional outputs are available from the ARKLS interface,
as listed in the following table and elaborated below.  We note that
where the name of an output would otherwise conflict with the
name of an optional output from the main solver, a suffix LS (for
Linear Solver) or MLS (for Mass Linear Solver) has been added here
(e.g. *lenrwLS*).



.. cssclass:: table-bordered

=================================================================  ========================================
Optional output                                                    Function name
=================================================================  ========================================
Stored Jacobian of the ODE RHS function                            :c:func:`ARKodeGetJac`
Time at which the Jacobian was evaluated                           :c:func:`ARKodeGetJacTime`
Step number at which the Jacobian was evaluated                    :c:func:`ARKodeGetJacNumSteps`
Size of real and integer workspaces                                :c:func:`ARKodeGetLinWorkSpace`
No. of Jacobian evaluations                                        :c:func:`ARKodeGetNumJacEvals`
No. of preconditioner evaluations                                  :c:func:`ARKodeGetNumPrecEvals`
No. of preconditioner solves                                       :c:func:`ARKodeGetNumPrecSolves`
No. of linear iterations                                           :c:func:`ARKodeGetNumLinIters`
No. of linear convergence failures                                 :c:func:`ARKodeGetNumLinConvFails`
No. of Jacobian-vector setup evaluations                           :c:func:`ARKodeGetNumJTSetupEvals`
No. of Jacobian-vector product evaluations                         :c:func:`ARKodeGetNumJtimesEvals`
No. of *fi* calls for finite diff. :math:`J` or :math:`Jv` evals.  :c:func:`ARKodeGetNumLinRhsEvals`
Last return from a linear solver function                          :c:func:`ARKodeGetLastLinFlag`
Name of constant associated with a return flag                     :c:func:`ARKodeGetLinReturnFlagName`
Size of real and integer mass matrix solver workspaces             :c:func:`ARKodeGetMassWorkSpace`
No. of mass matrix solver setups (incl. :math:`M` evals.)          :c:func:`ARKodeGetNumMassSetups`
No. of mass matrix multiply setups                                 :c:func:`ARKodeGetNumMassMultSetups`
No. of mass matrix multiplies                                      :c:func:`ARKodeGetNumMassMult`
No. of mass matrix solves                                          :c:func:`ARKodeGetNumMassSolves`
No. of mass matrix preconditioner evaluations                      :c:func:`ARKodeGetNumMassPrecEvals`
No. of mass matrix preconditioner solves                           :c:func:`ARKodeGetNumMassPrecSolves`
No. of mass matrix linear iterations                               :c:func:`ARKodeGetNumMassIters`
No. of mass matrix solver convergence failures                     :c:func:`ARKodeGetNumMassConvFails`
No. of mass-matrix-vector setup evaluations                        :c:func:`ARKodeGetNumMTSetups`
Last return from a mass matrix solver function                     :c:func:`ARKodeGetLastMassFlag`
=================================================================  ========================================

.. c:function:: int ARKodeGetJac(void* arkode_mem, SUNMatrix* J)

   Returns the internally stored copy of the Jacobian matrix of the ODE
   implicit right-hand side function.

   :param arkode_mem: the ARKODE memory structure.
   :param J: the Jacobian matrix.

   :retval ARKLS_SUCCESS: the output value has been successfully set.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

   .. warning::

      This function is provided for debugging purposes and the values in the
      returned matrix should not be altered.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetJacTime(void* arkode_mem, sunrealtype* t_J)

   Returns the time at which the internally stored copy of the Jacobian matrix
   of the ODE implicit right-hand side function was evaluated.

   :param arkode_mem: the ARKODE memory structure.
   :param t_J: the time at which the Jacobian was evaluated.

   :retval ARKLS_SUCCESS: the output value has been successfully set.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.


.. c:function:: int ARKodeGetJacNumSteps(void* arkode_mem, long int* nst_J)

   Returns the value of the internal step counter at which the internally stored copy of the
   Jacobian matrix of the ODE implicit right-hand side function was evaluated.

   :param arkode_mem: the ARKODE memory structure.
   :param nst_J: the value of the internal step counter at which the Jacobian was evaluated.

   :retval ARKLS_SUCCESS: the output value has been successfully set.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver interface has not been initialized.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetLinWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the real and integer workspace used by the ARKLS linear solver interface.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lenrwLS: the number of ``sunrealtype`` values in the ARKLS workspace.
   :param leniwLS: the number of integer values in the ARKLS workspace.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The workspace requirements reported by this routine
      correspond only to memory allocated within this interface and to
      memory allocated by the ``SUNLinearSolver`` object attached
      to it.  The template Jacobian matrix allocated by the user outside
      of ARKLS is not included in this report.

      In a parallel setting, the above values are global (i.e. summed over all
      processors).

   .. versionadded:: 6.1.0

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.


.. c:function:: int ARKodeGetNumJacEvals(void* arkode_mem, long int* njevals)

   Returns the number of Jacobian evaluations.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param njevals: number of Jacobian evaluations.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumPrecEvals(void* arkode_mem, long int* npevals)

   Returns the total number of preconditioner evaluations,
   i.e. the number of calls made to *psetup* with ``jok`` = ``SUNFALSE`` and
   that returned ``*jcurPtr`` = ``SUNTRUE``.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param npevals: the current number of calls to *psetup*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumPrecSolves(void* arkode_mem, long int* npsolves)

   Returns the number of calls made to the preconditioner
   solve function, *psolve*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param npsolves: the number of calls to *psolve*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumLinIters(void* arkode_mem, long int* nliters)

   Returns the cumulative number of linear iterations.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nliters: the current number of linear iterations.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumLinConvFails(void* arkode_mem, long int* nlcfails)

   Returns the cumulative number of linear convergence failures.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nlcfails: the current number of linear convergence failures.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumJTSetupEvals(void* arkode_mem, long int* njtsetup)

   Returns the cumulative number of calls made to the user-supplied
   Jacobian-vector setup function, *jtsetup*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param njtsetup: the current number of calls to *jtsetup*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumJtimesEvals(void* arkode_mem, long int* njvevals)

   Returns the cumulative number of calls made to the
   Jacobian-vector product function, *jtimes*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param njvevals: the current number of calls to *jtimes*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS)

   Returns the number of calls to the user-supplied implicit
   right-hand side function :math:`f^I` for finite difference
   Jacobian or Jacobian-vector product approximation.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nfevalsLS: the number of calls to the user implicit
                     right-hand side function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      The value *nfevalsLS* is incremented only if the default
      internal difference quotient function is used.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new linear solver
      module is "attached" to ARKODE, or when ARKODE is resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetLastLinFlag(void* arkode_mem, long int* lsflag)

   Returns the last return value from an ARKLS routine.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lsflag: the value of the last return flag from an
                  ARKLS function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: linear solvers are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.

      If the ARKLS setup function failed when using the
      ``SUNLINSOL_DENSE`` or ``SUNLINSOL_BAND`` modules, then the value
      of *lsflag* is equal to the column index (numbered from one) at
      which a zero diagonal element was encountered during the LU
      factorization of the (dense or banded) Jacobian matrix.  For all
      other failures, *lsflag* is negative.

      Otherwise, if the ARKLS setup function failed
      (:c:func:`ARKodeEvolve` returned *ARK_LSETUP_FAIL*), then
      *lsflag* will be *SUNLS_PSET_FAIL_UNREC*, *SUNLS_ASET_FAIL_UNREC*
      or *SUN_ERR_EXT_FAIL*.

      If the ARKLS solve function failed (:c:func:`ARKodeEvolve`
      returned *ARK_LSOLVE_FAIL*), then *lsflag* contains the error
      return flag from the ``SUNLinearSolver`` object, which will
      be one of:
      *SUN_ERR_ARG_CORRUPTRRUPT*, indicating that the ``SUNLinearSolver``
      memory is ``NULL``;
      *SUNLS_ATIMES_NULL*, indicating that a matrix-free iterative solver
      was provided, but is missing a routine for the matrix-vector product
      approximation,
      *SUNLS_ATIMES_FAIL_UNREC*, indicating an unrecoverable failure in
      the :math:`Jv` function;
      *SUNLS_PSOLVE_NULL*, indicating that an iterative linear solver was
      configured to use preconditioning, but no preconditioner solve
      routine was provided,
      *SUNLS_PSOLVE_FAIL_UNREC*, indicating that the preconditioner solve
      function failed unrecoverably;
      *SUNLS_GS_FAIL*, indicating a failure in the Gram-Schmidt procedure
      (SPGMR and SPFGMR only);
      *SUNLS_QRSOL_FAIL*, indicating that the matrix :math:`R` was found
      to be singular during the QR solve phase (SPGMR and SPFGMR only); or
      *SUN_ERR_EXT_FAIL*, indicating an unrecoverable failure in
      an external iterative linear solver package.

   .. versionadded:: 6.1.0


.. c:function:: char* ARKodeGetLinReturnFlagName(long int lsflag)

   Returns the name of the ARKLS constant corresponding to *lsflag*.

   :param lsflag: a return flag from an ARKLS function.

   :returns: The return value is a string containing the name of
             the corresponding constant. If using the ``SUNLINSOL_DENSE`` or
             ``SUNLINSOL_BAND`` modules, then if  1 :math:`\le` `lsflag`
             :math:`\le n` (LU factorization failed), this routine returns "NONE".

   .. warning::

      The user is responsible for freeing the returned string.

   .. note::

      This is only compatible with time-stepping modules that support implicit algebraic solvers.


   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS, long int* leniwMLS)

   Returns the real and integer workspace used by the ARKLS mass matrix linear solver interface.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lenrwMLS: the number of ``sunrealtype`` values in the ARKLS mass solver workspace.
   :param leniwMLS: the number of integer values in the ARKLS mass solver workspace.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      The workspace requirements reported by this routine
      correspond only to memory allocated within this interface and to
      memory allocated by the ``SUNLinearSolver`` object attached
      to it.  The template mass matrix allocated by the user outside
      of ARKLS is not included in this report.

      In a parallel setting, the above values are global (i.e. summed over all
      processors).

   .. versionadded:: 6.1.0

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.


.. c:function:: int ARKodeGetNumMassSetups(void* arkode_mem, long int* nmsetups)

   Returns the number of calls made to the ARKLS mass matrix solver
   'setup' routine; these include all calls to the user-supplied
   mass-matrix constructor function.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmsetups: number of calls to the mass matrix solver setup routine.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassMultSetups(void* arkode_mem, long int* nmvsetups)

   Returns the number of calls made to the ARKLS mass matrix 'matvec setup'
   (matrix-based solvers) routine.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmvsetups: number of calls to the mass matrix matrix-times-vector setup routine.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassMult(void* arkode_mem, long int* nmmults)

   Returns the number of calls made to the ARKLS mass matrix 'matvec'
   routine (matrix-based solvers) or the user-supplied *mtimes*
   routine (matris-free solvers).

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmmults: number of calls to the mass matrix solver matrix-times-vector routine.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassSolves(void* arkode_mem, long int* nmsolves)

   Returns the number of calls made to the ARKLS mass matrix solver 'solve' routine.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmsolves: number of calls to the mass matrix solver solve routine.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassPrecEvals(void* arkode_mem, long int* nmpevals)

   Returns the total number of mass matrix preconditioner evaluations,
   i.e. the number of calls made to *psetup*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmpevals: the current number of calls to *psetup*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassPrecSolves(void* arkode_mem, long int* nmpsolves)

   Returns the number of calls made to the mass matrix preconditioner
   solve function, *psolve*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmpsolves: the number of calls to *psolve*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassIters(void* arkode_mem, long int* nmiters)

   Returns the cumulative number of mass matrix solver iterations.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmiters: the current number of mass matrix solver linear iterations.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMassConvFails(void* arkode_mem, long int* nmcfails)

   Returns the cumulative number of mass matrix solver convergence failures.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmcfails: the current number of mass matrix solver convergence failures.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumMTSetups(void* arkode_mem, long int* nmtsetup)

   Returns the cumulative number of calls made to the user-supplied
   mass-matrix-vector product setup function, *mtsetup*.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nmtsetup: the current number of calls to *mtsetup*.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      This is only accumulated for the "life" of the linear
      solver object; the counter is reset whenever a new mass-matrix
      linear solver module is "attached" to ARKODE, or when ARKODE is
      resized.

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetLastMassFlag(void* arkode_mem, long int* mlsflag)

   Returns the last return value from an ARKLS mass matrix interface routine.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mlsflag: the value of the last return flag from an ARKLS
                   mass matrix solver interface function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARK_STEPPER_UNSUPPORTED: non-identity mass matrices are not supported
                                    by the current time-stepping module.

   .. note::

      This is only compatible with time-stepping modules that support non-identity mass matrices.

      The values of *msflag* for each of the various solvers
      will match those described above for the function
      :c:func:`ARKodeGetLastLinFlag`.

   .. versionadded:: 6.1.0




.. _ARKODE.Usage.ARKodeExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional routine may be called by a user to inquire
about existing solver parameters. While this would not typically be
called during the course of solving an initial value problem, it may
be useful for users wishing to better understand ARKODE.


.. cssclass:: table-bordered

====================================  ===================================
Optional routine                      Function name
====================================  ===================================
Output all ARKODE solver parameters   :c:func:`ARKodeWriteParameters`
====================================  ===================================




.. c:function:: int ARKodeWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all ARKODE solver parameters to the provided file pointer.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param fp: pointer to use for printing the solver parameters.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``.

   .. note::

      The *fp* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since parameters for all processes would be
      identical.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.write_parameters".

   .. versionadded:: 6.1.0


.. _ARKODE.Usage.Preallocation:

ARKODE data preallocation function
----------------------------------

Since the multi-stage structure of most ARKODE methods results in data
requirements that depend on the number of stages, ARKODE generally defers
allocation of stage-related internal data until the first call to
:c:func:`ARKodeEvolve`.  However, in some cases the user may wish to
preallocate this data earlier, for example to measure the memory footprint
before beginning a calculation, or to check for allocation errors at an
earlier time.  To request that that ARKODE preallocate all stage-related
internal data before the first call to :c:func:`ARKodeEvolve`, the user
may call the function :c:func:`ARKodeAllocateInternalData`.


.. c:function:: int ARKodeAllocateInternalData(void* arkode_mem)

   Optionally allocates stage-related internal data for the current ARKODE time-stepper module.

   :param arkode_mem: pointer to the ARKODE memory block.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_MEM_FAIL:  a memory allocation failed.

   .. warning::

      This must be called **after** all other optional input routines have been called,
      and **before** the first call to :c:func:`ARKodeEvolve`.  This routine should
      be called at most once per ARKODE memory block.

   .. versionadded:: x.y.z


.. _ARKODE.Usage.Reset:

ARKODE reset function
----------------------

To reset the ARKODE module to a particular state :math:`(t_R,y(t_R))` for the
continued solution of a problem, where a prior
call to ``*StepCreate`` has been made, the user must call the function
:c:func:`ARKodeReset`.  Like the stepper-specific ``*StepReInit`` functions,
this routine retains the current settings for all solver options and
performs no memory allocations but, unlike ``*StepReInit``, this routine
performs only a *subset* of the input checking and initializations that are
done in ``*StepCreate``. In particular this routine retains all internal
counter values and the step size/error history and does not reinitialize the
linear and/or nonlinear solver but it does indicate that a linear solver setup
is necessary in the next step. Like ``*StepReInit``, a call to
:c:func:`ARKodeReset` will delete any previously-set *tstop* value specified
via a call to :c:func:`ARKodeSetStopTime`.  Following a successful call to
:c:func:`ARKodeReset`, call :c:func:`ARKodeEvolve` again to continue
solving the problem. By default the next call to :c:func:`ARKodeEvolve` will
use the step size computed by ARKODE prior to calling :c:func:`ARKodeReset`.
To set a different step size or have ARKODE estimate a new step size use
:c:func:`ARKodeSetInitStep`.

One important use of the :c:func:`ARKodeReset` function is in the
treating of jump discontinuities in the RHS functions.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to :c:func:`ARKodeReset`.  To stop when
the location of the discontinuity is known, simply make that location
a value of ``tout``.  To stop when the location of the discontinuity
is determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS functions *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
functions (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int ARKodeReset(void* arkode_mem, sunrealtype tR, N_Vector yR)

   Resets the current ARKODE time-stepper module state to the provided
   independent variable value and dependent variable vector.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param tR: the value of the independent variable :math:`t`.
   :param yR: the value of the dependent variable vector :math:`y(t_R)`.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_MEM_FAIL:  a memory allocation failed.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      By default the next call to :c:func:`ARKodeEvolve` will use the step size
      computed by ARKODE prior to calling :c:func:`ARKodeReset`. To set a
      different step size or have ARKODE estimate a new step size use
      :c:func:`ARKodeSetInitStep`.

      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`ARKodeReset` also sends an error message to
      the error handler function.

   .. warning::

      Calling :c:func:`ARKodeReset` during forward integration of an IVP with
      checkpointing for adjoint sensitivity analysis is not supported.

   .. versionadded:: 6.1.0



.. _ARKODE.Usage.Resizing:

ARKODE system resize function
-------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatially-adaptive
PDE simulations under a method-of-lines approach), the ARKODE
integrator may be "resized" between integration steps, through calls
to the :c:func:`ARKodeResize` function. This function modifies
ARKODE's internal memory structures to use the new problem size,
without destruction of the temporal adaptivity heuristics.  It is
assumed that the dynamical time scales before and after the vector
resize will be comparable, so that all time-stepping heuristics prior
to calling :c:func:`ARKodeResize` remain valid after the call.  If
instead the dynamics should be recomputed from scratch, the ARKODE
memory structure should be deleted with a call to
:c:func:`ARKodeFree`, and recreated with a calls to
``*StepCreate``.

To aid in the vector resize operation, the user can supply a vector
resize function that will take as input a vector with the previous
size, and transform it in-place to return a corresponding vector of
the new size.  If this function (of type :c:func:`ARKVecResizeFn`)
is not supplied (i.e., is set to ``NULL``), then all existing vectors
internal to ARKODE will be destroyed and re-cloned from the new input
vector.

In the case that the dynamical time scale should be modified slightly
from the previous time scale, an input *hscale* is allowed, that will
rescale the upcoming time step by the specified factor.  If a value
*hscale* :math:`\le 0` is specified, the default of 1.0 will be used.



.. c:function:: int ARKodeResize(void* arkode_mem, N_Vector yR, sunrealtype hscale, sunrealtype tR, ARKVecResizeFn resize, void* resize_data)

   Re-sizes ARKODE with a different state vector but with comparable
   dynamical time scale.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param yR: the newly-sized state vector, holding the current
              dependent variable values :math:`y(t_R)`.
   :param hscale: the desired time step scaling factor (i.e. the next
                  step will be of size *h\*hscale*).
   :param tR: the current value of the independent variable
              :math:`t_R` (this must be consistent with *yR*).
   :param resize: the user-supplied vector resize function (of type
                  :c:func:`ARKVecResizeFn`.
   :param resize_data: the user-supplied data structure to be passed
                       to *resize* when modifying internal ARKODE vectors.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_NULL:  ``arkode_mem`` was ``NULL``.
   :retval ARK_NO_MALLOC: ``arkode_mem`` was not allocated.
   :retval ARK_ILL_INPUT: an argument had an illegal value.

   .. note::

      If an error occurred, :c:func:`ARKodeResize` also sends an error
      message to the error handler function.

      If inequality constraint checking is enabled a call to
      :c:func:`ARKodeResize` will disable constraint checking. A call
      to :c:func:`ARKodeSetConstraints` is required to re-enable constraint
      checking.

      **Resizing the linear solver:**

      When using any of the SUNDIALS-provided linear solver modules, the
      linear solver memory structures must also be resized.  At present,
      none of these include a solver-specific "resize" function, so the linear
      solver memory must be destroyed and re-allocated **following** each
      call to :c:func:`ARKodeResize`.  Moreover, the existing ARKLS
      interface should then be deleted and recreated by attaching the
      updated ``SUNLinearSolver`` (and possibly ``SUNMatrix``) object(s)
      through calls to
      :c:func:`ARKodeSetLinearSolver`, and
      :c:func:`ARKodeSetMassLinearSolver`.

      If any user-supplied routines are provided to aid the linear solver
      (e.g. Jacobian construction, Jacobian-vector product,
      mass-matrix-vector product, preconditioning), then the corresponding
      "set" routines must be called again **following** the solver
      re-specification.

      **Resizing the absolute tolerance array:**

      If using array-valued absolute tolerances, the absolute tolerance
      vector will be invalid after the call to :c:func:`ARKodeResize`, so
      the new absolute tolerance vector should be re-set **following** each
      call to :c:func:`ARKodeResize` through a new call to
      :c:func:`ARKodeSVtolerances` and possibly
      :c:func:`ARKodeResVtolerance` if applicable.

      If scalar-valued tolerances or a tolerance function was specified
      through either :c:func:`ARKodeSStolerances` or
      :c:func:`ARKodeWFtolerances`, then these will remain valid and no
      further action is necessary.

      **Example codes:**

      * ``examples/arkode/C_serial/ark_heat1D_adapt.c``

   .. versionadded:: 6.1.0


.. _ARKODE.Usage.MRIStepInterface:

Using an ARKODE solver as an MRIStep "inner" solver
---------------------------------------------------

When using an integrator from ARKODE as the inner (fast) integrator with MRIStep, the
utility function :c:func:`ARKodeCreateMRIStepInnerStepper` should be used to
wrap the ARKODE memory block as an :c:type:`MRIStepInnerStepper`.

.. c:function:: int ARKodeCreateMRIStepInnerStepper(void *inner_arkode_mem, MRIStepInnerStepper *stepper)

   Wraps an ARKODE integrator as an :c:type:`MRIStepInnerStepper` for use
   with MRIStep.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param stepper: the :c:type:`MRIStepInnerStepper` object to create.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_FAIL: a memory allocation failed.
   :retval ARK_STEPPER_UNSUPPORTED: the time-stepping module does not currently support use as an inner stepper.

   .. note::

      Currently, ARKODE integrators based on ARKStep, ERKStep, and MRIStep
      support use as an MRIStep inner stepper.

   **Example usage:**

      .. code-block:: C

         /* fast (inner) and slow (outer) ARKODE objects */
         void *inner_arkode_mem = NULL;
         void *outer_arkode_mem = NULL;

         /* MRIStepInnerStepper to wrap the inner (fast) object */
         MRIStepInnerStepper stepper = NULL;

         /* create an ARKODE object, setting fast (inner) right-hand side
            functions and the initial condition */
         inner_arkode_mem = *StepCreate(...);

         /* configure the inner integrator */
         retval = ARKodeSet*(inner_arkode_mem, ...);

         /* create MRIStepInnerStepper wrapper for the ARKODE integrator */
         flag = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

         /* create an MRIStep object, setting the slow (outer) right-hand side
            functions and the initial condition */
         outer_arkode_mem = MRIStepCreate(fse, fsi, t0, y0, stepper, sunctx)


.. _ARKODE.Usage.SUNStepperInterface:

Using an ARKODE solver as a SUNStepper
--------------------------------------

The utility function :c:func:`ARKodeCreateSUNStepper` wraps an ARKODE memory
block as a :c:type:`SUNStepper`.

.. c:function:: int ARKodeCreateSUNStepper(void *inner_arkode_mem, SUNStepper *stepper)

   Wraps an ARKODE integrator as a :c:type:`SUNStepper`.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param stepper: the :c:type:`SUNStepper` object.

   :retval ARK_SUCCESS: the function exited successfully.
   :retval ARK_MEM_FAIL: a memory allocation failed.
   :retval ARK_SUNSTEPPER_ERR: the :c:type:`SUNStepper` initialization failed.

   .. warning::
      Currently, ``stepper`` will be equipped with an implementation for the
      :c:func:`SUNStepper_SetForcing` function only if ``inner_arkode_mem`` is
      an ARKStep, ERKStep, or MRIStep integrator.

   .. versionadded:: 6.2.0
