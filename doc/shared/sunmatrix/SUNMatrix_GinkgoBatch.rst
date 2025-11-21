.. _SUNMatrix.GinkgoBatch:

The SUNMATRIX_GINKGOBATCH Module
================================

.. versionadded:: 7.5.0

The SUNMATRIX_GINKGOBATCH implementation of the ``SUNMatrix`` API provides an
interface to the batched matrix types from the Ginkgo linear algebra library.
This module is written in C++17 and is distributed as a header file. To use the
SUNMATRIX_GINKGOBATCH ``SUNMatrix``, users will need to include
``sunmatrix/sunmatrix_ginkgobatch.hpp``. The module is meant to be used with the
SUNLINEARSOLVER_GINKGOBATCH module described in :numref:`SUNLinSol.GinkgoBatch`.

.. note::

   It is assumed that users of this module are aware of how to use Ginkgo. This
   module does not try to encapsulate Ginkgo matrices, rather it provides a
   lightweight interoperability layer between Ginkgo and SUNDIALS. Most, if not
   all, of the Ginkgo batch matrix types should work with this interface.

.. _SUNMatrix.GinkgoBatch.CompatibleNVectors:

Compatible Vectors
------------------

The :c:type:`N_Vector` to use with the SUNLINEARSOLVER_GINKGOBATCH module depends on the ``gko::Executor``
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

Compatible Packages
-------------------

This module will work with any of the SUNDIALS packages. The only caveat is that, when using ARKODE with a
non-identity mass matrix, the only Ginkgo matrix type currently supported is ``BatchDense``.


.. _SUNMatrix.GinkgoBatch.API:

SUNMATRIX_GINKGOBATCH API
-------------------------

In this section we list the public API of the
:cpp:type:`sundials::ginkgo::BatchMatrix` class.

.. cpp:class:: template<class GkoBatchMatType> sundials::ginkgo::BatchMatrix : public sundials::ConvertibleTo<SUNMatrix>

   Batched matrix wrapper for Ginkgo batch matrix types, providing a SUNDIALS
   ``SUNMatrix`` interface.

   .. cpp:function:: BatchMatrix()

      Default constructor. The matrix must be copied or moved to.

   .. cpp:function:: BatchMatrix(gko::size_type num_batches, sunindextype M, sunindextype N, std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)

      Construct a batch matrix with the given number of batches, rows ``M``, columns ``N``,
      executor, and context. (Specialized for supported Ginkgo batch matrix
      types.)

   .. cpp:function:: BatchMatrix(gko::size_type num_batches, sunindextype M, sunindextype N, sunindextype num_nonzeros, std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)

      Construct a batch sparse matrix with the given number of batches, rows ``M``,
      columns ``N``, nonzeros, executor, and context. (Specialized for supported
      Ginkgo batch matrix types.)

   .. cpp:function:: BatchMatrix(std::shared_ptr<GkoBatchMatType> gko_mat, SUNContext sunctx)

      Construct a BatchMatrix from an existing Ginkgo batch matrix pointer and
      SUNDIALS context.

   .. cpp:function:: BatchMatrix(BatchMatrix&& that_matrix) noexcept

      Move constructor.

   .. cpp:function:: BatchMatrix(const BatchMatrix& that_matrix)

      Copy constructor. Clones the Ginkgo matrix and SUNDIALS SUNMatrix.

   .. cpp:function:: BatchMatrix& operator=(BatchMatrix&& rhs) noexcept

      Move assignment.

   .. cpp:function:: BatchMatrix& operator=(const BatchMatrix& rhs)

      Copy assignment. Clones the Ginkgo matrix and SUNDIALS SUNMatrix.

   .. cpp:function:: ~BatchMatrix() override = default

      Default destructor.

   .. cpp:function:: std::shared_ptr<GkoBatchMatType> GkoMtx() const

      Get the underlying Ginkgo batch matrix pointer.

   .. cpp:function:: std::shared_ptr<const gko::Executor> GkoExec() const

      Get the Ginkgo executor associated with the matrix.

   .. cpp:function:: const gko::batch_dim<2>& GkoSize() const

      Get the Ginkgo batch size object.

   .. cpp:function:: sunindextype NumBatches() const

      Get the number of batches (batch systems).

   .. cpp:function:: operator SUNMatrix() override

      Implicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: operator SUNMatrix() const override

      Implicit conversion to a :c:type:`SUNMatrix`.

   .. cpp:function:: SUNMatrix get() override

      Explicit conversion to a :c:type:`SUNMatrix`.

      .. versionadded:: x.y.z
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: SUNMatrix get() const override

      Explicit conversion to a :c:type:`SUNMatrix`.

      .. versionadded:: x.y.z
      
         Replaces the ``Convert`` method which was deprecated.
