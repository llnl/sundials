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

.. _SUNMatrix.Ginkgo:

The SUNMATRIX_GINKGO Module
===========================

.. versionadded:: 6.4.0

The SUNMATRIX_GINKGO implementation of the ``SUNMatrix`` API provides an
interface to the matrix data structure for the Ginkgo linear algebra library :cite:p:`ginkgo-toms-2022`.
Ginkgo provides several different matrix formats and linear solvers which can run on a variety
of hardware, such as NVIDIA, AMD, and Intel GPUs as well as multicore CPUs.
Since Ginkgo is a modern C++ library, SUNMATRIX_GINKGO is also written in
modern C++ (it requires C++14). Unlike most other SUNDIALS modules, it is a header only library.
To use the SUNMATRIX_GINKGO ``SUNMatrix``, users will need to include ``sunmatrix/sunmatrix_ginkgo.hpp``.
More instructions on building SUNDIALS with Ginkgo enabled are given in
:numref:`Installation.Options.Ginkgo`. For instructions on building and using
Ginkgo itself, refer to the `Ginkgo website and documentation <https://ginkgo-project.github.io/>`_.

.. note::

  It is assumed that users of this module are aware of how to use Ginkgo. This module does not
  try to encapsulate Ginkgo matrices, rather it provides a lightweight iteroperability layer
  between Ginkgo and SUNDIALS.

The SUNMATRIX_GINKGO module is defined by the ``sundials::ginkgo::Matrix`` templated class:

.. code-block:: cpp

  template<typename GkoMatType>
  class Matrix : public sundials::impl::BaseMatrix, public sundials::ConvertibleTo<SUNMatrix>;

.. _SUNMatrix.Ginkgo.CompatibleNVectors:

Compatible Vectors
------------------

The :c:type:`N_Vector` to use with the SUNLINEARSOLVER_GINKGO module depends on the ``gko::Executor``
utilized. That is, when using the ``gko::CudaExecutor`` you should use a CUDA capable ``N_Vector``
(e.g., :numref:`NVectors.CUDA`), ``gko::HipExecutor`` goes with a HIP capable ``N_Vector`` (e.g.,
:numref:`NVectors.HIP`), ``gko::DpcppExecutor`` goes with a DPC++/SYCL capable ``N_Vector`` (e.g.,
:numref:`NVectors.SYCL`),  and ``gko::OmpExecutor`` goes with a CPU based N_Vector (e.g.,
:numref:`NVectors.OpenMP`). Specifically, what makes a ``N_Vector`` compatible with different Ginkgo
executors is where they store the data. The GPU enabled Ginkgo executors need the data to reside on
the GPU, so the ``N_Vector`` must implement :c:func:`N_VGetDeviceArrayPointer` and keep the data in
GPU memory. The CPU-only enabled Ginkgo executors (e.g, ``gko::OmpExecutor`` and
``gko::ReferenceExecutor``) need data to reside on the CPU and will use
:c:func:`N_VGetArrayPointer` to access the ``N_Vector`` data.

.. _SUNMatrix.Ginkgo.Usage:

Using SUNMATRIX_GINKGO
----------------------

To use the SUNMATRIX_GINKGO module, we begin by creating an instance of a Ginkgo matrix using
Ginkgo's API. For example, below we create a Ginkgo sparse matrix that uses the CSR storage format
and then fill the diagonal of the matrix with ones to make an identity matrix:

.. code-block:: cpp

   auto gko_matrix{gko::matrix::Csr<sunrealtype, sunindextype>::create(gko_exec, matrix_dim)};
   gko_matrix->read(gko::matrix_data<sunrealtype, sunindextype>::diag(matrix_dim, 1.0));

After we have a Ginkgo matrix object, we wrap it in an instance of the ``sundials::ginkgo::Matrix``
class. This object can be provided to other SUNDIALS functions that expect a ``SUNMatrix`` object
via implicit conversion, or the ``get()`` method:

.. code-block:: cpp

  sundials::ginkgo::Matrix<gko::matrix::Csr> matrix{gko_matrix, sunctx};
  SUNMatrix I1 = matrix.get(); // explicit conversion to SUNMatrix
  SUNMatrix I2 = matrix;  // implicit conversion to SUNMatrix

No further interaction with ``matrix`` is required from this point, and it is possible to
to use the ``SUNMatrix`` API operating on ``I1`` or ``I2`` (or if needed, via Ginkgo operations
on ``gko_matrix``).


.. warning::

  :c:func:`SUNMatDestroy` should never be called on a ``SUNMatrix`` that was created via conversion
  from a ``sundials::ginkgo::Matrix``. Doing so may result in a double free.


.. _SUNMatrix.Ginkgo.API:

SUNMATRIX_GINKGO API
--------------------

In this section we list the public API of the ``sundials::ginkgo::Matrix`` class.

.. cpp:class:: template<typename GkoMatType> \
               Matrix : public sundials::impl::BaseMatrix, public sundials::ConvertibleTo<SUNMatrix>

  .. cpp:function:: Matrix() = default

      Default constructor - means the matrix must be copied or moved to.

  .. cpp:function:: Matrix(std::shared_ptr<GkoMatType> gko_mat, SUNContext sunctx)

      Constructs a Matrix from an existing Ginkgo matrix object.

      :param gko_mat:  A Ginkgo matrix object
      :param sunctx: The SUNDIALS simulation context object (:c:type:`SUNContext`)

  .. cpp:function:: Matrix(Matrix&& that_matrix) noexcept

      Move constructor.

  .. cpp:function:: Matrix(const Matrix& that_matrix)

      Copy constructor (performs a deep copy).

  .. cpp:function:: Matrix& operator=(Matrix&& rhs) noexcept

      Move assignment.

  .. cpp:function:: Matrix& operator=(const Matrix& rhs)

      Copy assignment clones the ``gko::matrix`` and :c:type:`SUNMatrix`.
      This is a deep copy (i.e. a new data array is created).

  .. cpp:function:: virtual ~Matrix() = default;

      Default destructor.

  .. cpp:function:: std::shared_ptr<GkoMatType> GkoMtx() const

      Get the underlying Ginkgo matrix object.

  .. cpp:function:: std::shared_ptr<const gko::Executor> GkoExec() const

      Get the ``gko::Executor`` associated with the Ginkgo matrix.

  .. cpp:function:: const gko::dim<2>& GkoSize() const

      Get the size, i.e. ``gko::dim``, for the Ginkgo matrix.

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
