..
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

.. _SUNLinSol.GinkgoBatch:

The SUNLINEARSOLVER_GINKGOBATCH Module
=======================================

.. versionadded:: 7.5.0

The SUNLINEARSOLVER_GINKGOBATCH implementation of the ``SUNLinearSolver`` API provides an
interface to the batched linear solvers from the Ginkgo linear algebra library :cite:p:`ginkgo-toms-2022`.
Like SUNLINEARSOLVER_GINKGO, this module is written in C++17 and is distributed as a header file.
To use the SUNLINEARSOLVER_GINKGOBATCH ``SUNLinearSolver``, users will need to include 
``sunlinsol/sunlinsol_ginkgobatch.hpp``. The module is meant to be used with the SUNMATRIX_GINKGOBATCH 
module described in :numref:`SUNMatrix.GinkgoBatch`. Instructions on building SUNDIALS  with Ginkgo enabled are given
in :numref:`Installation.Options.Ginkgo`.  For instructions on building and using Ginkgo itself, refer to the
`Ginkgo website and documentation <https://ginkgo-project.github.io/>`__.

.. note::

   It is assumed that users of this module are aware of how to use Ginkgo. This module does not
   try to encapsulate Ginkgo linear solvers, rather it provides a lightweight iteroperability layer
   between Ginkgo and SUNDIALS. Most, if not all, of the Ginkgo linear solvers should work with this
   interface.

.. _SUNLinSol.GinkgoBatch.Usage:

Using SUNLINEARSOLVER_GINKGOBATCH
---------------------------------

After choosing a compatible ``N_Vector`` (see :numref:`SUNMatrix.GinkgoBatch.CompatibleNVectors`) and creating a Ginkgo-enabled ``SUNMatrix`` (see
:numref:`SUNMatrix.GinkgoBatch`) to use the SUNLINEARSOLVER_GINKGOBATCH module, we create the linear solver object:

.. code-block:: cpp

   using GkoBatchMatrixType = gko::batch::matrix::Csr<sunrealtype, sunindextype>;
   using GkoBatchSolverType = gko::batch::solver::Bicgstab<sunrealtype>;
   using SUNGkoMatrixType   = sundials::ginkgo::BatchMatrix<GkoBatchMatrixType>;
   using SUNGkoLinearSolverType =
      sundials::ginkgo::BatchLinearSolver<GkoBatchSolverType, GkoBatchMatrixType>;

   SUNGkoLinearSolverType LS{gko_exec, gko::batch::stop::tolerance_type::absolute,
                             precond_factory, num_batches, sunctx};

Next, we can pass the instance of ``sundials::ginkgo::BatchLinearSolver`` to any function
expecting a ``SUNLinearSolver`` object through the implicit conversion operator or explicit conversion function.
For example,

.. code-block:: cpp

   // Attach linear solver and matrix to CVODE.
   //
   // Implicit conversion from sundials::ginkgo::BatchLinearSolver<GkoBatchSolverType, GkoBatchMatrixType>
   // to a SUNLinearSolver object is done.
   //
   // For details about creating A see the SUNMATRIX_GINKGOBATCH module.
   CVodeSetLinearSolver(cvode_mem, LS, A);

   // Alternatively with explicit conversion of LS to a SUNLinearSolver
   // and A to a SUNMatrix:
   CVodeSetLinearSolver(cvode_mem, LS.get(), A.get());

After attaching the linear solver to the SUNDIALS integrator, one must change the norm factor the integrator uses
since the Ginkgo linear solver will take norms over individual batches, not the entire system.

.. code-block:: cpp

   // When using ARKODE:
   ARKodeSetLSNormFactor(arkode_mem, std::sqrt(batch_size));

   // When using CVODE:
   CVodeSetLSNormFactor(cvode_mem, std::sqrt(batch_size));

   // When using IDA:
   IDASetLSNormFactor(ida_mem, std::sqrt(batch_size));

.. warning:: 

   Setting the linear solver norm factor is essential. If this is not set, you will likely see a large 
   number of linear solver convergence failures.

.. warning::

   :c:func:`SUNLinSolFree` should never be called on a ``SUNLinearSolver`` that was created via conversion
   from a :cpp:type:`sundials::ginkgo::BatchLinearSolver`. Doing so may result in a double free.


.. _SUNLinSol.GinkgoBatch.API:

SUNLINEARSOLVER_GINKGOBATCH API
-------------------------------

All `core functions <SUNLinSol.CoreFn>` of the ``SUNLinearSolver`` API are supported by this module/class.
However, we note a difference in behavior for :c:func:`SUNLinSolNumIters`:

.. c:function:: int SUNLinSolNumIters_GinkgoBatch(SUNLinearSolver S)

   This function returns the average number of iterations across all of the
   batch systems. As such, functions that utilize this function to accumulate
   statistics over steps or solve (i.e., the ``GetNumLinIters`` function in each
   package) will return the sum of these average values.


The public API of the :cpp:type:`sundials::ginkgo::BatchLinearSolver` class is as follows:

.. cpp:class:: template<class GkoBatchSolverType, class GkoBatchMatType> \
               sundials::ginkgo::BatchLinearSolver : public sundials::ConvertibleTo<SUNLinearSolver>

   .. cpp::member:: NO_SCALING

   .. cpp::member:: LAGGED_SCALING

   .. cpp::member:: SOLVE_SCALING

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with default tolerance type and max iterations.

      :param gko_exec: The `gko::Executor` to use
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with specified tolerance type.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with a preconditioner factory.

      :param gko_exec: The `gko::Executor` to use
      :param precon_factory: Ginkgo batch preconditioner factory
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, int max_iters, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with a maximum number of iterations.

      :param gko_exec: The `gko::Executor` to use
      :param max_iters: Maximum number of iterations
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, int max_iters, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with specified tolerance type and maximum iterations.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param max_iters: Maximum number of iterations
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, int max_iters, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with a preconditioner factory and maximum iterations.

      :param gko_exec: The `gko::Executor` to use
      :param precon_factory: Ginkgo batch preconditioner factory
      :param max_iters: Maximum number of iterations
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with specified tolerance type and preconditioner factory.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param precon_factory: Ginkgo batch preconditioner factory
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::batch::stop::tolerance_type tolerance_type, std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory, int max_iters, sunindextype num_batches, SUNContext sunctx)

      Constructs a new BatchLinearSolver with all options specified.

      :param gko_exec: The `gko::Executor` to use
      :param tolerance_type: Ginkgo batch solver tolerance type
      :param precon_factory: Ginkgo batch preconditioner factory
      :param max_iters: Maximum number of iterations
      :param num_batches: Number of batches (batch systems)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: BatchLinearSolver(BatchLinearSolver&& that_solver) noexcept

      Move constructor.

   .. cpp:function:: BatchLinearSolver& operator=(BatchLinearSolver&& rhs)

      Move assignment.

   .. cpp:function:: ~BatchLinearSolver() override = default

      Default destructor.

   .. cpp:function:: operator SUNLinearSolver() override

      Implicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: operator SUNLinearSolver() const override

      Implicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: SUNLinearSolver get() override

      Explicit conversion to a :c:type:`SUNLinearSolver`.

      .. versionadded:: 7.6.0
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: SUNLinearSolver get() const override

      Explicit conversion to a :c:type:`SUNLinearSolver`.

      .. versionadded:: 7.6.0
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: std::shared_ptr<const gko::Executor> GkoExec() const

      Get the ``gko::Executor`` associated with the Ginkgo solver.

   .. cpp:function:: std::shared_ptr<typename GkoBatchSolverType::Factory> GkoFactory()

      Get the underlying Ginkgo solver factory.

   .. cpp:function:: GkoBatchSolverType* GkoSolver()

      Get the underlying Ginkgo solver.

      .. note::

         This will be ``nullptr`` until the linear solver setup phase.

   .. cpp:function:: int AvgNumIters() const

      Get the average number of linear solver iterations across the batches in the most recent solve.

   .. cpp:function:: int StddevNumIters() const

      Get the standard deviation of the number of iterations across the batches during the last solve.

   .. cpp:function:: int SumAvgNumIters() const

      Get the running sum of the average number of iterations in this solver's lifetime.

   .. cpp:function:: SUNErrCode SetScalingMode(int scaling_mode)

      Sets the matrix scaling mode. The options are:

      * ``BatchLinearSolver::NO_SCALING`` -- no scaling of the matrix
      * ``BatchLinearSolver::LAGGED_SCALING`` -- the matrix is only scaled when it is updated, this is the default
      * ``BatchLinearSolver::SOLVE_SCALING`` -- the matrix is scaled (and unscaled) every solve, this is the most expensive option on a per-solve basis

   .. cpp:function:: SUNErrCode SetScalingVectors(N_Vector s1, N_Vector s2)

      Sets the left (``s1``) and right (``s2``) scaling vectors to use.

   .. cpp:function:: int Setup(BatchMatrix<GkoBatchMatType>* A)

      Setup the linear system.

      :param A: the linear system matrix

      :returns: Zero on success otherwise a non-zero value

   .. cpp:function:: int Solve(BatchMatrix<GkoBatchMatType>* A, N_Vector b, N_Vector x, sunrealtype tol)

      Solve the linear system Ax = b to the specified tolerance.

      :param A: the linear system matrix
      :param b: the right-hand side vector
      :param x: the solution vector
      :param tol: the tolerance to solve the system to

      :returns: Zero on success otherwise a non-zero value
