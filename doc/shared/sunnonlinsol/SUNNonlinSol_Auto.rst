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

============================================
The SUNNonlinSol_Auto implementation
============================================

This section describes the SUNNonlinSol implementation that can automatically
switch between :numref:`SUNNonlinSol.FixedPoint` and :numref:`SUNNonlinSol.Newton`
during a solve.

To access the SUNNonlinSol_Auto module, include the header file
``sunnonlinsol/sunnonlinsol_auto.h``. The library to link to is
``libsundials_sunnonlinsolauto.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.


.. _SUNNonlinSol.Auto.Description:

SUNNonlinSol_Auto description
----------------------------------------

SUNNonlinSol_Auto is a hybrid nonlinear solver that delegates each nonlinear
solve to an underlying fixed-point or Newton solver and may request a switch
to the other algorithm based on a runtime switching criterion.

Switching decisions are checked at every nonlinear iteration by intercepting
the convergence test function provided by the integrator (see
:c:func:`SUNNonlinSolSetConvTestFn`). When SUNNonlinSol_Auto decides to switch,
it returns :c:macro:`SUN_NLS_SWITCH` from the convergence test. The underlying
solver propagates this code back to the integrator. Integrators that support
automatic nonlinear solver switching will respond by reinitializing the
nonlinear solver interface and retrying the nonlinear solve.


.. _SUNNonlinSol.Auto.Functions:

SUNNonlinSol_Auto functions
----------------------------------------

The SUNNonlinSol_Auto module provides the following constructor for creating the
``SUNNonlinearSolver`` object.

.. c:function:: SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m, SUNNonlinSolAutoType active_solver_type, SUNContext sunctx)

   This creates a ``SUNNonlinearSolver`` object for use with SUNDIALS integrators.

   **Arguments:**
      * *y* -- a template for cloning vectors needed within the solver.
      * *m* -- the number of acceleration vectors to use with the underlying
        fixed-point solver (passed to :c:func:`SUNNonlinSol_FixedPoint`).
      * *active_solver_type* -- the initial solver active_solver_type (see :c:active_solver_type:`SUNNonlinSolAutoType`).
      * *sunctx* -- the :c:active_solver_type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      A SUNNonlinSol object if the constructor exits successfully, otherwise it
      will be ``NULL``.

.. c:enum:: SUNNonlinSolAutoType

   An identifier indicating which underlying solver is used initially.

   .. c:enumerator:: SUNNONLINSOL_AUTO_FIXEDPOINT

      Use the fixed-point solver.

   .. c:enumerator:: SUNNONLINSOL_AUTO_NEWTON

      Use the Newton solver.

The SUNNonlinSol_Auto module implements all of the functions defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn` except for
:c:func:`SUNNonlinSolSetup`. The SUNNonlinSol_Auto functions have the same names
as those defined by the generic SUNNonlinSol API with ``_Auto`` appended to the
function name.

The SUNNonlinSol_Auto module also defines the following user-callable functions
to configure the switching criteria:

.. c:function:: SUNErrCode SUNNonlinSolSetFpToNewtAlpha_Auto(SUNNonlinearSolver NLS, sunrealtype alpha)

   Set the fixed-point criterion parameter ``alpha`` where the switching test is
   ``crate >= alpha`` (with ``crate`` computed by the fixed-point solver).

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *alpha* -- the criterion parameter (must be positive).

   **Return value:**
      A :c:active_solver_type:`SUNErrCode`.

.. c:function:: SUNErrCode SUNNonlinSolSetNewtToFpThreshold_Auto(SUNNonlinearSolver NLS, sunrealtype threshold)

   Set the Newton criterion parameter ``threshold`` where the switching test is
   ``stiffr < threshold`` (with ``stiffr`` computed by the Newton solver).

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *threshold* -- the criterion parameter (must be positive).

   **Return value:**
      A :c:active_solver_type:`SUNErrCode`.

