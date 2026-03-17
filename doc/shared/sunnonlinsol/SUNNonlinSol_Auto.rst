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
it returns :c:macro:`SUN_NLS_SWITCH` from the convergence test. The underlying
solver propagates this code back to the integrator. Integrators that support
automatic nonlinear solver switching will respond by reinitializing the
nonlinear solver interface and retrying the nonlinear solve.

A full mathematical description of the switching criterion and algorithm can be
found in :cite:p:`norsett1986switching`. In short, switching from fixed-point to 
Newton occurs when the rate of convergence of the solver, defined as

.. math:: 

   R_ \leftarrow \max\{0.3R, \|\delta_m\| / \|\delta_{m-1}\|\},

is indicating slow convergence or divergence, i.e., when ``R > \alpha``, where 
:math:`0.0 < \alpha < 1.0`, for a specified number of consecutive iterations (default is 1).
Switching from Newton to fixed-point occurs when the stiffness indicator, defined as 

.. math::

   \text{stiffr} \leftarrow \|F(y^n)\| / \|\delta_{m} \|,

satisfies :math:`\text{stiffr} < \beta` where :math:`1.0 < \beta \leq 2` for a certain number
of consecutive iterations (default is 10).



.. _SUNNonlinSol.Auto.Functions:

SUNNonlinSol_Auto functions
---------------------------

The SUNNonlinSol_Auto module provides the following constructor for creating the
``SUNNonlinearSolver`` object.

.. c:function:: SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m, SUNNonlinSolAutoType active_solver_type, SUNContext sunctx)

   This creates a ``SUNNonlinearSolver`` object for use with SUNDIALS integrators.

   :param y: a template for cloning vectors needed within the solver.
   :param m: the number of acceleration vectors to use with the underlying fixed-point solver (passed to :c:func:`SUNNonlinSol_FixedPoint`).
   :param active_solver_type: the initial solver active_solver_type (see :c:active_solver_type:`SUNNonlinSolAutoType`).
   :param sunctx: the :c:active_solver_type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :returns: a pointer to a SUNNonlinSol object if the constructor exits successfully, otherwise it will be ``NULL``.

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

The SUNNonlinSol_Auto module also defines the following functions for controlling the switching behavior.

.. c:function:: SUNErrCode SUNNonlinSolSetSwitchingParameters_Auto(SUNNonlinearSolver NLS, sunrealtype newt_to_fp_threshold, \
                                                                   long int newt_to_fp_delay, sunrealtype fp_to_newt_threshold, \
                                                                   long int fp_to_newt_delay)

   This function sets the parameters that control the switching behavior of the SUNNonlinSol_Auto module. 

   :param NLS: the nonlinear solver object returned by :c:func:`SUNNonlinSol_Auto`.
   :param newt_to_fp_threshold: the threshold for switching from Newton to fixed-point (aka :math:`\beta`)
   :param newt_to_fp_delay: the number of consecutive iterations that must satisfy the switching criterion 
      before switching from Newton to fixed-point.
   :param fp_to_newt_threshold: the threshold for switching from fixed-point to Newton (aka :math:`\alpha`).
   :param fp_to_newt_delay: the number of consecutive iterations that must satisfy the switching 
      criterion before switching from fixed-point to Newton.

   :returns: :c:macro:`SUN_SUCCESS` if successful, otherwise an error code.
