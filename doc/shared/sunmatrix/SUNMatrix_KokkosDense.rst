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

.. _SUNMatrix.Kokkos:

The SUNMATRIX_KOKKOSDENSE Module
================================

.. versionadded:: 6.4.0

The SUNMATRIX_KOKKOSDENSE :c:type:`SUNMatrix` implementation provides a data
structure for dense and dense batched (block-diagonal) matrices using Kokkos
:cite:p:`edwards2014kokkos,trott2022kokkos` and KokkosKernels
:cite:p:`trott2021kokkos` to support a variety of backends including serial, OpenMP,
CUDA, HIP, and SYCL. Since Kokkos is a modern C++ library, the module is also
written in modern C++ (it requires C++14) as a header only library. To utilize
this ``SUNMatrix`` users will need to include
``sunmatrix/sunmatrix_kokkosdense.hpp``. More instructions on building SUNDIALS
with Kokkos and KokkosKernels enabled are given in
:numref:`Installation.Options.KokkosKernels`. For instructions on building and
using Kokkos and KokkosKernels, refer to the
`Kokkos <https://kokkos.github.io/kokkos-core-wiki/index.html>`_
and `KokkosKernels <https://github.com/kokkos/kokkos-kernels/wiki>`_.
documentation.

.. _SUNMatrix.Kokkos.Usage:

Using SUNMATRIX_KOKKOSDENSE
----------------------------

The SUNMATRIX_KOKKOSDENSE module is defined by the ``DenseMatrix`` templated
class in the ``sundials::kokkos`` namespace:

.. code-block:: cpp

   template<class ExecutionSpace = Kokkos::DefaultExecutionSpace,
            class MemorySpace = typename ExecutionSpace::memory_space>
   class DenseMatrix : public sundials::impl::BaseMatrix,
                       public sundials::ConvertibleTo<SUNMatrix>

To use the SUNMATRIX_KOKKOSDENSE module, we begin by constructing an instance of
the Kokkos dense matrix e.g.,

.. code-block:: cpp

   // Single matrix using the default execution space
   sundials::kokkos::DenseMatrix<> A{rows, cols, sunctx};

   // Batched (block-diagonal) matrix using the default execution space
   sundials::kokkos::DenseMatrix<> Abatch{blocks, rows, cols, sunctx};

   // Batched (block-diagonal) matrix using the Cuda execution space
   sundials::kokkos::DenseMatrix<Kokkos::Cuda> Abatch{blocks, rows, cols, sunctx};

   // Batched (block-diagonal) matrix using the Cuda execution space and
   // a non-default execution space instance
   sundials::kokkos::DenseMatrix<Kokkos::Cuda> Abatch{blocks, rows, cols,
                                                      exec_space_instance,
                                                      sunctx};

Instances of the ``DenseMatrix`` class are implicitly or explicitly (using the
:cpp:func:`~DenseMatrix::get` method) convertible to a :c:type:`SUNMatrix`
e.g.,

.. code-block:: cpp

   sundials::kokkos::DenseMatrix<> A{rows, cols, sunctx};
   SUNMatrix B = A;           // implicit conversion to SUNMatrix
   SUNMatrix C = A.get(); // explicit conversion to SUNMatrix

No further interaction with a ``DenseMatrix`` is required from this point, and
it is possible to use the :c:type:`SUNMatrix` API to operate on ``B`` or ``C``.

.. warning::

   :c:func:`SUNMatDestroy` should never be called on a ``SUNMatrix`` that was
   created via conversion from a ``sundials::kokkos::DenseMatrix``. Doing so may
   result in a double free.

The underlying ``DenseMatrix`` can be extracted from a ``SUNMatrix`` using
:cpp:func:`GetDenseMat` e.g.,

.. code-block:: cpp

   auto A_dense_mat = GetDenseMat<>(A_sunmat);

The SUNMATRIX_KOKKOSDENSE module is compatible with the NVECTOR_KOKKOS vector
module (see :numref:`NVectors.Kokkos`) and SUNLINEARSOLVER_KOKKOSDENSE linear
solver module (see :numref:`SUNLinSol.Kokkos`).


.. _SUNMatrix.Kokkos.API:

SUNMATRIX_KOKKOSDENSE API
-------------------------

In this section we list the public API of the ``sundials::kokkos::DenseMatrix``
class.

.. cpp:class:: template<class ExeccutionSpace = Kokkos::DefaultExecutionSpace, \
                        class MemorySpace = typename ExecutionSpace::memory_space> \
               DenseMatrix : public sundials::impl::BaseMatrix, \
                             public sundials::ConvertibleTo<SUNMatrix>

   .. cpp:type:: exec_space   = ExecutionSpace;
   .. cpp:type:: memory_space = MemorySpace;
   .. cpp:type:: view_type    = Kokkos::View<sunrealtype***, memory_space>;
   .. cpp:type:: size_type    = typename view_type::size_type;
   .. cpp:type:: range_policy = Kokkos::MDRangePolicy<exec_space, Kokkos::Rank<3>>;
   .. cpp:type:: team_policy  = typename Kokkos::TeamPolicy<exec_space>;
   .. cpp:type:: member_type  = typename Kokkos::TeamPolicy<exec_space>::member_type;

   .. cpp:function:: DenseMatrix() = default

      Default constructor -- the matrix must be copied or moved to.

   .. cpp:function:: DenseMatrix(size_type rows, size_type cols, \
                                 SUNContext sunctx)

      Constructs a single DenseMatrix using the default execution space
      instance.

      :param rows: number of matrix rows
      :param cols: number of matrix columns
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(size_type rows, size_type cols, \
                                 exec_space ex, SUNContext sunctx)

      Constructs a single DenseMatrix using the provided execution space
      instance.

      :param rows: number of matrix rows
      :param cols: number of matrix columns
      :param ex: an execution space
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(size_type blocks, size_type block_rows, \
                                 size_type block_cols, SUNContext sunctx)

      Constructs a batched (block-diagonal) DenseMatrix using the default
      execution space instance.

      :param blocks: number of matrix blocks
      :param block_rows: number of rows in a block
      :param block_cols: number of columns in a block
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(size_type blocks, size_type block_rows, \
                                 size_type block_cols, exec_space ex, \
                                 SUNContext sunctx)

      Constructs a batched (block-diagonal) DenseMatrix using the provided
      execution space instance.

      :param blocks: number of matrix blocks
      :param block_rows: number of rows in a block
      :param block_cols: number of columns in a block
      :param ex: an execution space
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: DenseMatrix(DenseMatrix&& that_matrix) noexcept

      Move constructor.

   .. cpp:function:: DenseMatrix(const DenseMatrix& that_matrix)

      Copy constructor. This creates a shallow clone of the Matrix, i.e., it
      creates a new Matrix with the same properties, such as size, but it does
      not copy the data.

   .. cpp:function:: DenseMatrix& operator=(DenseMatrix&& rhs) noexcept

      Move assignment.

   .. cpp:function:: DenseMatrix& operator=(const DenseMatrix& rhs)

      Copy assignment. This creates a shallow clone of the Matrix, i.e., it
      creates a new Matrix with the same properties, such as size, but it does
      not copy the data.

   .. cpp:function:: virtual ~DenseMatrix() = default;

      Default destructor.

   .. cpp:function:: exec_space ExecSpace()

      Get the execution space instance used by the matrix.

   .. cpp:function:: view_type View()

      Get the underlying Kokkos view with extents
      ``{blocks, block_rows, block_cols}``.

   .. cpp:function:: size_type Blocks()

      Get the number of blocks i.e., ``extent(0)``.

   .. cpp:function:: size_type BlockRows()

      Get the number of rows in a block i.e., ``extent(1)``.

   .. cpp:function:: size_type BlockCols()

      Get the number of columns in a block i.e., ``extent(2)``.

   .. cpp:function:: size_type Rows()

      Get the number of rows in the block-diagonal matrix i.e.,
      ``extent(0) * extent(1)``.

   .. cpp:function:: size_type Cols()

      Get the number of columns in the block-diagonal matrix i.e.,
      ``extent(0) * extent(2)``.

   .. cpp:function:: operator SUNMatrix() override

      Implicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: operator SUNMatrix() const override

      Implicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: SUNMatrix get() override

      Explicit conversion to a :c:type:`SUNMatrix`.

      .. versionadded:: 7.6.0
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: SUNMatrix get() const override

      Explicit conversion to a :c:type:`SUNMatrix`.

      .. versionadded:: 7.6.0
      
         Replaces the ``Convert`` method which was deprecated.

.. cpp:function:: template<class ExecutionSpace = Kokkos::DefaultExecutionSpace, \
                           class MemorySpace = typename ExecutionSpace::memory_space> \
                  inline DenseMatrix<MatrixType>* GetDenseMat(SUNMatrix A)

   Get the dense matrix wrapped by a SUNMatrix
