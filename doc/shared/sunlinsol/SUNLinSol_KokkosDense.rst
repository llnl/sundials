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

.. _SUNLinSol.Kokkos:

The SUNLINEARSOLVER_KOKKOSDENSE Module
======================================

.. versionadded:: 6.4.0

The SUNLINEARSOLVER_KOKKOSDENSE :c:type:`SUNLinearSolver` implementation
provides an interface to KokkosKernels :cite:p:`trott2021kokkos` linear solvers
for dense and batched dense (block-diagonal) systems. Since Kokkos is a modern
C++ library, the module is also written in modern C++ (it requires C++14) as a
header only library. To utilize this ``SUNLinearSolver`` user will need to
include ``sunlinsol/sunlinsol_kokkosdense.hpp``. More instructions on building
SUNDIALS with Kokkos and KokkosKernels enabled are given in
:numref:`Installation.Options.KokkosKernels`. For instructions on building and
using Kokkos and KokkosKernels, refer to the
`Kokkos <https://kokkos.github.io/kokkos-core-wiki/index.html>`_
and `KokkosKernels <https://github.com/kokkos/kokkos-kernels/wiki>`_.
documentation.

.. _SUNLinSol.Kokkos.Usage:

Using SUNLINEARSOLVER_KOKKOSDENSE
---------------------------------

The SUNLINEARSOLVER_KOKKOSDENSE module is defined by the ``DenseLinearSolver``
templated class in the ``sundials::kokkos`` namespace:

.. code-block:: cpp

   template<class ExecSpace = Kokkos::DefaultExecutionSpace,
            class MemSpace = typename ExecSpace::memory_space>
   class DenseLinearSolver : public sundials::impl::BaseLinearSolver,
                             public sundials::ConvertibleTo<SUNLinearSolver>

To use the SUNLINEARSOLVER_KOKKOSDENSE module, we begin by constructing an
instance of a dense linear solver e.g.,

.. code-block:: cpp

    // Create a dense linear solver
    sundials::kokkos::DenseLinearSolver<> LS{sunctx};

Instances of the ``DenseLinearSolver`` class are implicitly or explicitly (using
the :cpp:func:`~DenseLinearSolver::get` method) convertible to a
:c:type:`SUNLinearSolver` e.g.,

.. code-block:: cpp

   sundials::kokkos::DenseLinearSolver<> LS{sunctx};
   SUNLinearSolver LSA = LS;           // implicit conversion to SUNLinearSolver
   SUNLinearSolver LSB = LS.get(); // explicit conversion to SUNLinearSolver

.. warning::

  :c:func:`SUNLinSolFree` should never be called on a ``SUNLinearSolver`` that
  was created via conversion from a ``sundials::kokkos::DenseLinearSolver``.
  Doing so may result in a double free.

The SUNLINEARSOLVER_KOKKOSDENSE module is compatible with the NVECTOR_KOKKOS
vector module (see :numref:`NVectors.Kokkos`) and SUNMATRIX_KOKKOSDENSE matrix
module (see :numref:`SUNMatrix.Kokkos`).


.. _SUNLinSol.Kokkos.API:

SUNLINEARSOLVER_KOKKOSDENSE API
-------------------------------

In this section we list the public API of the
``sundials::kokkos::DenseLinearSolver`` class.

.. cpp:class:: template<class ExecSpace = Kokkos::DefaultExecutionSpace, \
                        class MemSpace = typename ExecSpace::memory_space> \
               DenseLinearSolver : public sundials::impl::BaseLinearSolver, \
                                   public sundials::ConvertibleTo<SUNLinearSolver>

   .. cpp:function:: DenseLinearSolver() = default;

      Default constructor - means the solver must be moved to.

   .. cpp:function:: DenseLinearSolver(SUNContext sunctx)

      Constructs a new DenseLinearSolver.

      :param sunctx: The SUNDIALS simulation context (:c:type:`SUNContext`)

   .. cpp:function:: DenseLinearSolver(DenseLinearSolver&& that_solver) noexcept

      Move constructor.

   .. cpp:function:: DenseLinearSolver& operator=(DenseLinearSolver&& rhs)

      Move assignment.

   .. cpp:function:: ~DenseLinearSolver() override = default

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
