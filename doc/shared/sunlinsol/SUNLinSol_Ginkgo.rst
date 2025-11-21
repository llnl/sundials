..
   ----------------------------------------------------------------
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
   ----------------------------------------------------------------

.. _SUNLinSol.Ginkgo:

The SUNLINEARSOLVER_GINKGO Module
=================================

.. versionadded:: 6.4.0

The SUNLINEARSOLVER_GINKGO implementation of the ``SUNLinearSolver`` API provides an
interface to the linear solvers from the Ginkgo linear algebra library :cite:p:`ginkgo-toms-2022`.
Since Ginkgo is a modern C++ library, SUNLINEARSOLVER_GINKGO is also written in
modern C++ (specifically, C++14). Unlike most other SUNDIALS modules, it is
a header only library. To use the SUNLINEARSOLVER_GINKGO ``SUNLinearSolver``, users will
need to include ``sunlinsol/sunlinsol_ginkgo.hpp``. The module is meant to be used with
the SUNMATRIX_GINKGO module described in :numref:`SUNMatrix.Ginkgo`.
Instructions on building SUNDIALS  with Ginkgo enabled are given
in :numref:`Installation.Options.Ginkgo`.  For instructions on
building and using Ginkgo itself, refer to the
`Ginkgo website and documentation <https://ginkgo-project.github.io/>`_.

.. note::

  It is assumed that users of this module are aware of how to use Ginkgo. This module does not
  try to encapsulate Ginkgo linear solvers, rather it provides a lightweight iteroperability layer
  between Ginkgo and SUNDIALS. Most, if not all, of the Ginkgo linear solver should work with this
  interface.

.. _SUNLinSol.Ginkgo.Usage:

Using SUNLINEARSOLVER_GINKGO
----------------------------

After choosing a compatible ``N_Vector`` (see :numref:`SUNMatrix.Ginkgo.CompatibleNVectors`) and creating a Ginkgo-enabled ``SUNMatrix`` (see
:numref:`SUNMatrix.Ginkgo`) to use the SUNLINEARSOLVER_GINKGO module, we first create a Ginkgo
stopping criteria object. Importantly, the ``sundials::ginkgo::DefaultStop`` class provided
by SUNDIALS implements a stopping criterion that matches the default SUNDIALS stopping criterion.
Namely, it checks if the max iterations (5 by default) were reached or if the absolute residual
norm was below a specified tolerance. The criterion can be created just like any other
Ginkgo stopping criteria:

.. code-block:: cpp

   auto crit{sundials::ginkgo::DefaultStop::build().with_max_iters(max_iters).on(gko_exec)};

.. warning::
   It is *highly* recommended to employ this criterion when using Ginkgo solvers with SUNDIALS,
   but it is optional. However, to use the Ginkgo multigrid or cbgmres linear solvers, different
   Ginkgo criterion must be used.

Once we have created our stopping criterion, we create a Ginkgo solver factory object and
wrap it in a :cpp:type:`sundials::ginkgo::LinearSolver` object. In this example, we create
a Ginkgo conjugate gradient solver:

.. code-block:: cpp

   using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
   using GkoSolverType = gko::solver::Cg<sunrealtype>;

   auto gko_solver_factory = gko::share(
      GkoSolverType::build().with_criteria(std::move(crit)).on(gko_exec));

   sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType> LS{
      gko_solver_factory, sunctx};

Finally, we can pass the instance of :cpp:type:`sundials::ginkgo::LinearSolver` to any function
expecting a ``SUNLinearSolver`` object through the implicit conversion operator or explicit conversion function.

.. code-block:: cpp

   // Attach linear solver and matrix to CVODE.
   //
   // Implicit conversion from sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>
   // to a SUNLinearSolver object is done.
   //
   // For details about creating A see the SUNMATRIX_GINKGO module.
   CVodeSetLinearSolver(cvode_mem, LS, A);

   // Alternatively with explicit conversion of LS to a SUNLinearSolver
   // and A to a SUNMatrix:
   CVodeSetLinearSolver(cvode_mem, LS->get(), A->get());


.. warning::

  :c:func:`SUNLinSolFree` should never be called on a ``SUNLinearSolver`` that was created via conversion
  from a :cpp:type:`sundials::ginkgo::LinearSolver`. Doing so may result in a double free.


.. _SUNLinSol.Ginkgo.API:

SUNLINEARSOLVER_GINKGO API
--------------------------

In this section we list the public API of the :cpp:type:`sundials::ginkgo::LinearSolver` class.

.. cpp:class:: template<class GkoSolverType, class GkoMatrixType> \
               sundials::ginkgo::LinearSolver : public sundials::ConvertibleTo<SUNLinearSolver>

   .. cpp:function:: LinearSolver() = default;

      Default constructor - means the solver must be moved to.

   .. cpp:function:: LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory, SUNContext sunctx)

      Constructs a new LinearSolver from a Ginkgo solver factory.

      :param gko_solver_factory: The Ginkgo solver factory (typically `gko::matrix::<type>::Factory``)
      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: LinearSolver(LinearSolver&& that_solver) noexcept

      Move constructor.

   .. cpp:function:: LinearSolver& operator=(LinearSolver&& rhs)

      Move assignment.

   .. cpp:function:: ~LinearSolver() override = default

      Default destructor.

   .. cpp:function:: operator SUNLinearSolver() override

      Implicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: operator SUNLinearSolver() const override

      Implicit conversion to a :c:type:`SUNLinearSolver`.

   .. cpp:function:: SUNLinearSolver get() override

      Explicit conversion to a :c:type:`SUNLinearSolver`.

      .. versionadded:: x.y.z
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: SUNLinearSolver get() const override

      Explicit conversion to a :c:type:`SUNLinearSolver`.

      .. versionadded:: x.y.z
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: std::shared_ptr<const gko::Executor> GkoExec() const

      Get the ``gko::Executor`` associated with the Ginkgo solver.

   .. cpp:function:: std::shared_ptr<typename GkoSolverType::Factory> GkoFactory()

      Get the underlying Ginkgo solver factory.

   .. cpp:function:: GkoSolverType* GkoSolver()

      Get the underlying Ginkgo solver.

      .. note::

         This will be `nullptr` until the linear solver setup phase.

   .. cpp:function:: int NumIters() const

      Get the number of linear solver iterations in the most recent solve.

   .. cpp:function:: sunrealtype ResNorm() const

      Get the residual norm of the solution at the end of the last solve.

      The type of residual norm depends on the Ginkgo stopping criteria
      used with the solver. With the ``DefaultStop`` criteria this would
      be the absolute residual 2-norm.

   .. cpp:function:: GkoSolverType* Setup(Matrix<GkoMatrixType>* A)

      Setup the linear system.

      :param A: the linear system matrix

      :returns: Pointer to the Ginkgo solver generated from the factory

   .. cpp:function:: gko::LinOp* Solve(N_Vector b, N_Vector x, sunrealtype tol)

      Solve the linear system Ax = b to the specified tolerance.

      :param b: the right-hand side vector
      :param x: the solution vector
      :param tol: the tolerance to solve the system to

      :returns: ``gko::LinOp*`` the solution
