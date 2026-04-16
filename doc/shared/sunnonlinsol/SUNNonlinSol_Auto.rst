..
   Programmer(s): Cody J. Balos @ LLNL
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

.. _SUNNonlinSol.Auto:

====================================
The SUNNonlinSol_Auto implementation
====================================

This section describes the SUNNonlinSol implementation that can automatically
switch between :numref:`SUNNonlinSol.FixedPoint` and :numref:`SUNNonlinSol.Newton`
during a solve. The switching algorithm is based on :cite:p:`norsett1986switching`.

To access the SUNNonlinSol_Auto module, include the header file
``sunnonlinsol/sunnonlinsol_auto.h``. The library to link to is
``libsundials_sunnonlinsolauto.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.


.. _SUNNonlinSol.Auto.Description:

SUNNonlinSol_Auto description
-----------------------------

SUNNonlinSol_Auto is a hybrid nonlinear solver that delegates each nonlinear
solve to an underlying fixed-point or Newton solver and may request a switch
to the other algorithm based on a runtime switching criterion.

Switching decisions are checked at every nonlinear iteration by intercepting
the convergence test function provided by the integrator (see
:c:func:`SUNNonlinSolSetConvTestFn`). When SUNNonlinSol_Auto decides to switch,
it returns ``SUN_NLS_SWITCH`` from the convergence test. The underlying
solver propagates this code back to the integrator. Integrators that support
automatic nonlinear solver switching will respond by reinitializing the
nonlinear solver interface and retrying the nonlinear solve.

A full mathematical description of the switching criterion and algorithm can be
found in :cite:p:`norsett1986switching`. In short, switching from fixed-point to
Newton occurs when the fixed-point convergence-rate estimate

.. math::

   R \leftarrow \max\{0.3R, \|\delta_m\| / \|\delta_{m-1}\|\},

indicates slow convergence or divergence, i.e., when :math:`R \ge \alpha`, for
a specified number of nonlinear solves. The implementation default is
:math:`\alpha = 0.8` with no additional delay for switching from fixed-point to
Newton. Switching from Newton to fixed-point occurs when the stiffness
indicator

.. math::

   \text{stiffr} \leftarrow \|F(y^m)\| / \|\delta_m\|,

satisfies :math:`\text{stiffr} < \beta` for a specified number of nonlinear
solves. The implementation default is :math:`\beta = 2.0` with a delay of
10 solves before switching from Newton to fixed-point is allowed.



.. _SUNNonlinSol.Auto.Functions:

SUNNonlinSol_Auto functions
---------------------------

The SUNNonlinSol_Auto module provides the following constructor for creating the
``SUNNonlinearSolver`` object.

.. c:function:: SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m, SUNNonlinSolAutoType initial_solver_type, SUNContext sunctx)

   This creates a ``SUNNonlinearSolver`` object for use with SUNDIALS integrators.

   :param y: a template for cloning vectors needed within the solver.
   :param m: the number of acceleration vectors to use with the underlying fixed-point solver
      (passed to :c:func:`SUNNonlinSol_FixedPoint`).
   :param initial_solver_type: the initial solver type.
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :returns: a pointer to a ``SUNNonlinearSolver`` object if the constructor exits successfully, otherwise it will be ``NULL``.

.. c:enum:: SUNNonlinSolAutoType

   An identifier indicating which underlying solver is used initially.

   .. c:enumerator:: SUNNONLINSOL_AUTO_FIXEDPOINT

      Use the fixed-point solver.

   .. c:enumerator:: SUNNONLINSOL_AUTO_NEWTON

      Use the Newton solver.

The SUNNonlinSol_Auto module implements the generic nonlinear-solver interface
defined in :numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`
except for :c:func:`SUNNonlinSolSetup`. Users should normally call the generic
SUNNonlinSol API. When solver-specific access is needed, the auto module
provides ``_Auto`` helper routines for switching parameters, accumulated
statistics, and access to the underlying fixed-point and Newton solvers.

The SUNNonlinSol_Auto module also defines the following function for
controlling the switching behavior.

.. c:function:: SUNErrCode SUNNonlinSolSetSwitchingParameters_Auto(SUNNonlinearSolver NLS, sunrealtype newt_to_fp_threshold, long int newt_to_fp_delay, sunrealtype fp_to_newt_threshold, long int fp_to_newt_delay)

   This function sets the parameters that control the switching behavior of the
   ``SUNNonlinearSolver_Auto`` module.

   :param NLS: the nonlinear solver object returned by :c:func:`SUNNonlinSol_Auto`.
   :param newt_to_fp_threshold: the threshold for switching from Newton to fixed-point (i.e., :math:`\beta`).
   :param newt_to_fp_delay: the minimum number of nonlinear solves after a switch before switching from Newton to fixed-point is allowed.
   :param fp_to_newt_threshold: the threshold for switching from fixed-point to Newton (i.e., :math:`\alpha`).
   :param fp_to_newt_delay: the minimum number of nonlinear solves after a switch before switching from fixed-point to Newton is allowed.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an error code.

   .. note::

      If supported by the SUNNonlinearSolver implementation, this routine will be called by
      :c:func:`SUNNonlinSolSetOptions` when using the keys
      ``NLSid.newt_to_fp_threshold``, ``NLSid.newt_to_fp_delay``,
      ``NLSid.fp_to_newt_threshold``, and ``NLSid.fp_to_newt_delay``.

.. c:function:: SUNErrCode SUNNonlinSolGetFixedPointSolver_Auto(SUNNonlinearSolver NLS, SUNNonlinearSolver *fp_nls)

   This function returns the underlying fixed-point solver so that users may
   configure fixed-point-specific options directly, e.g.,
   :c:func:`SUNNonlinSolSetMaxIters` or
   :c:func:`SUNNonlinSolSetConvRateConstant_FixedPoint`.

   :param NLS: a ``SUNNonlinearSolver`` object returned by :c:func:`SUNNonlinSol_Auto`.
   :param fp_nls: a pointer to the underlying fixed-point solver object.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an error code.

.. c:function:: SUNErrCode SUNNonlinSolGetNewtonSolver_Auto(SUNNonlinearSolver NLS, SUNNonlinearSolver *newton_nls)

   This function returns the underlying Newton solver so that users may
   configure Newton-specific options directly, e.g.,
   :c:func:`SUNNonlinSolSetMaxIters`.

   :param NLS: a ``SUNNonlinearSolver`` object returned by :c:func:`SUNNonlinSol_Auto`.
   :param newton_nls: a pointer to the underlying Newton solver object.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an error code.

.. c:function:: SUNErrCode SUNNonlinSolGetTotalNumItersByType_Auto(SUNNonlinearSolver NLS, long int *fp_iters, long int *newt_iters)

   This function returns the total nonlinear iteration counts accumulated by the
   fixed-point and Newton sub-solvers.

   :param NLS: a ``SUNNonlinearSolver`` object returned by :c:func:`SUNNonlinSol_Auto`.
   :param fp_iters: the total number of fixed-point iterations.
   :param newt_iters: the total number of Newton iterations.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an error code.

.. c:function:: SUNErrCode SUNNonlinSolGetNumConvFailsByType_Auto(SUNNonlinearSolver NLS, long int *fp_nconvfails, long int *newt_nconvfails)

   This function returns the total nonlinear convergence failures accumulated by
   the fixed-point and Newton sub-solvers.

   :param NLS: a ``SUNNonlinearSolver`` object returned by :c:func:`SUNNonlinSol_Auto`.
   :param fp_nconvfails: the total number of fixed-point convergence failures.
   :param newt_nconvfails: the total number of Newton convergence failures.

   :returns: ``SUN_SUCCESS`` if successful, otherwise an error code.
