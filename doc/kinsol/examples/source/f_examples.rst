..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _KINSOL.Examples.Fortran:

Fortran Example Problems
========================

The Fortran example problem programs supplied with the KINSOL package are all
written using the Fortran 2003 standard, use double precision arithmetic, and
32- or 64-bit integers for indexing. See the :ref:`Fortran` section in
the KINSOL User Documentation for details.

A serial example: kinDiagon_kry_f2003
-------------------------------------

The ``kinDiagon_kry_f2003`` program solves the problem

.. math::

   F(u) = 0, \quad \text{where} \quad F_i(u) = u_i^2 - i^2, \quad \text{for} \quad 1 \le i \le N.

using the :ref:`serial vector <NVectors.NVSerial>` and the :ref:`GMRES linear
solver <SUNLinSol.SPGMR>` with a diagonal preconditioner.

The ``kinDiagonKry_mod`` module includes the problem parameters (e.g., ``neq``
:math:`=N`) used in subroutines and functions contained in the module. The
subroutines and function defined in the module are:

* The subroutine ``init`` sets the components of the initial guess vector
  (``sunvec_u``) to the values :math:`u_i = 2 i`, all the components of the
  scaling vector (``sunvec_s``) are set to ``1.0`` (i.e., no scaling is done),
  and the components of the constraints vector (``sunvec_c``) are set to ``0.0``
  (i.e., no inequality constraints are imposed on the solution vector).

* The function ``func`` defines the user-supplied nonlinear residual
  :math:`F(u)`.

* The functions ``kpsetup`` and ``kpsolve`` are the preconditioner setup and
  solve functions, respectively. ``kpsetup`` fills the array ``p`` with an
  approximation to the reciprocals of the Jacobian diagonal elements which are
  then used in ``kpsolve`` to solve the preconditioner linear system :math:`P x
  = v` by simple multiplication.

The ``main`` program starts with a number of ``use`` lines, which allow access
to KINSOL (``fkinsol_mod``), the serial vector (``fnvector_serial_mod``), and
the GMRES linear solver (``fsunlinsol_spgmr_mod``). This is followed by a number
of problem parameters for configuring the linear and nonlinear solver. In
particular, the maximum number of iterations between calls to the preconditioner
setup routine (``msbpre = 5``), the tolerance for stopping based on the function
norm (``fnormtol = 1.0d-5``), and the tolerance for stopping based on the step
length (``scsteptol = 1.0d-4``) are specified. After the problem parameters are
declarations for local variables including SUNDIALS vector, matrix, and linear
solver objects.

After printing the problem description, the ``main`` program creates the
simulation context (``sunctx``) that is passed to all SUNDIALS constructors.
The program creates then creates serial solution, scaling, and constraint
vectors using ``FN_VNew_Serial`` and sets the vector values by calling the
``init`` subroutine.

The KINSOL solver is created by calling ``FKINCreate`` and initialized with
``FKINInit`` which takes as input the ``c_funloc`` of the Fortran residual
function (``func``), a vector to use as a template to create internal workspace
(``sunvec_u``), and the SUNDIALS context. Then solver options are specified by
calling various ``FKINSet`` functions. In particular, the maximum number of
iterations between calls to the preconditioner setup routine
(``FKINSetMaxSetupCalls``), the tolerance for stopping based on the function
norm (``FKINSetFuncNormTol``), and the tolerance for stopping based on the step
length (``FKINSetScaledStepTol``).

Next, the GMRES linear solver is created with ``FSUNLinSol_SPGMR`` specifying
the maximum Krylov subspace dimension (``maxl = 10``) and right preconditioning
(``prectype = SUN_PREC_RIGHT``). The linear solver is then attached to KINSOL
with ``FKINSetLinearSolver``. The input matrix (``sunmat_J``) is assigned to
``null`` for the matrix-free linear solver. The maximum number of restarts
allowed for GMRES is then updated to ``maxlrst = 2`` by calling
``FSUNLinSol_SPGMRSetMaxRestarts``. The preconditioner functions are then
attached by calling ``FKINSetPreconditioner`` and passing the ``c_funloc`` of
``kpsetup`` and ``kpsolve``.

The solution of the nonlinear system is obtained after a successful return from
``FKINSol``, which is then printed using the ``PrintOutput`` subroutine. Solver
statistics are then output using the ``PrintFinalStats`` function which calls
various ``FKINGet`` functions. Finally, the memory allocated for the KINSOL, the
linear solver, vectors, and context are released by calling ``FKINFree``,
``FSUNLinSolFree``, ``FN_VDestroy``, and ``FSUNContext_Free``, respectively.

The following is sample output from ``kinDiagon_kry_f2003``.

.. literalinclude:: ../../../../examples/kinsol/F2003_serial/kinDiagon_kry_f2003.out
   :language: none

A parallel example: kin_diagon_kry_f2003
----------------------------------------

The program ``kin_diagon_kry_f2003`` is a straightforward modification of
``kinDiagon_kry_f2003`` to use the MPI vector.

After initialing MPI, the MPI parallel vectors are created using
``FN_VNew_Parallel`` with the default MPI communicator and the local and global
vector sizes as its first three arguments. The rank of the local process,
``myid``, is used in both the initial guess and the system function, inasmuch as
the global and local indices to the vector ``u`` are related by ``iglobal =
ilocal + mype*nlocal``. In other respects, the problem setup (KINSOL
initialization, linear solver creation and specification, etc.) and solution
steps are the same as in ``kinDiagon_kry_f2003``. Upon successful return from
``FKINSol``, the solution on all ranks is printed and statistics are output from
the root processes. Finally, the allocated memory is released and the MPI
environment is finalized.

For this simple example, no inter-process communication is required to evaluate
the nonlinear system function or the preconditioner. As a consequence, the
user-supplied routines ``func``, ``kpsetup``, and ``kpsolve`` are basically
identical to those in ``kinDiagon_kry_f2003``.

The following is sample output from ``kin_diagon_kry_2003`` using 4 MPI ranks.

.. literalinclude:: ../../../../examples/kinsol/F2003_parallel/kin_diagon_kry_f2003.out
   :language: none
