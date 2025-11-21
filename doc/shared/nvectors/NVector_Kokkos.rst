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

.. _NVectors.Kokkos:

The NVECTOR_KOKKOS Module
=========================

.. versionadded:: 6.4.0

The NVECTOR_KOKKOS :c:type:`N_Vector` implementation provides a vector data
structure using Kokkos :cite:p:`edwards2014kokkos,trott2022kokkos` to support a
variety of backends including serial, OpenMP, CUDA, HIP, and SYCL. Since Kokkos is
a modern C++ library, the module is also written in modern C++ (it requires
C++14) as a header only library. To utilize this ``N_Vector`` users will need to
include ``nvector/nvector_kokkos.hpp``. More instructions on building SUNDIALS
with Kokkos enabled are given in :numref:`Installation.Options.Kokkos`.
For instructions on building and using Kokkos, refer to the `Kokkos
<https://kokkos.github.io/kokkos-core-wiki/index.html>`_ documentation.

.. _NVectors.Kokkos.Usage:

Using NVECTOR_KOKKOS
--------------------

The NVECTOR_KOKKOS module is defined by the ``Vector`` templated
class in the ``sundials::kokkos`` namespace:

.. code-block:: cpp

   template<class ExecutionSpace = Kokkos::DefaultExecutionSpace,
            class MemorySpace = typename ExecutionSpace::memory_space>
      class Vector : public sundials::impl::BaseNVector,
                     public sundials::ConvertibleTo<N_Vector>

To use the NVECTOR_KOKKOS module, we construct an instance of the ``Vector`` class e.g.,

.. code-block:: cpp

   // Vector with extent length using the default execution space
   sundials::kokkos::Vector<> x{length, sunctx};

   // Vector with extent length using the Cuda execution space
   sundials::kokkos::Vector<Kokkos::Cuda> x{length, sunctx};

   // Vector based on an existing Kokkos::View
   Kokkos::View<> view{"a view", length};
   sundials::kokkos::Vector<> x{view, sunctx};

   // Vector based on an existing Kokkos::View for device and host
   Kokkos::View<Kokkos::Cuda> device_view{"a view", length};
   Kokkos::View<Kokkos::HostMirror> host_view{Kokkos::create_mirror_view(device_view)};
   sundials::kokkos::Vector<> x{device_view, host_view, sunctx};


Instances of the ``Vector`` class are implicitly or explicitly (using the
:cpp:func:`~Vector::get` method) convertible to a :c:type:`N_Vector`
e.g.,

.. code-block:: cpp

   sundials::kokkos::Vector<> x{length, sunctx};
   N_Vector x2 = x;       // implicit conversion to N_Vector
   N_Vector x3 = x.get(); // explicit conversion to N_Vector

No further interaction with a ``Vector`` is required from this point, and
it is possible to use the :c:type:`N_Vector` API to operate on ``x2`` or ``x3``.

.. warning::

   :c:func:`N_VDestroy` should never be called on a ``N_Vector`` that was
   created via conversion from a ``sundials::kokkos::Vector``. Doing so may
   result in a double free.

The underlying ``Vector`` can be extracted from a ``N_Vector`` using
:cpp:func:`GetVec` e.g.,

.. code-block:: cpp

   auto x_vec = GetVec<>(x3);

.. _NVectors.Kokkos.API:

NVECTOR_KOKKOS API
------------------

In this section we list the public API of the ``sundials::kokkos::Vector``
class.

.. cpp:class:: template<class ExecutionSpace = Kokkos::DefaultExecutionSpace, \
                        class MemorySpace = class ExecutionSpace::memory_space> \
               Vector : public sundials::impl::BaseNVector, \
                        public sundials::ConvertibleTo<N_Vector>

   .. cpp:type:: view_type      = Kokkos::View<sunrealtype*, MemorySpace>;
   .. cpp:type:: size_type      = typename view_type::size_type;
   .. cpp:type:: host_view_type = typename view_type::HostMirror;
   .. cpp:type:: memory_space   = MemorySpace;
   .. cpp:type:: exec_space     = typename MemorySpace::execution_space;
   .. cpp:type:: range_policy   = Kokkos::RangePolicy<exec_space>;

   .. cpp:function:: Vector() = default

      Default constructor -- the vector must be copied or moved to.

   .. cpp:function:: Vector(size_type length, SUNContext sunctx)

      Constructs a single ``Vector`` which is based on a 1D ``Kokkos::View``
      with the ExecutionSpace and MemorySpace provided as template arguments.

      :param length: length of the vector (i.e., the extent of the View)
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: Vector(view_type view, SUNContext sunctx)

      Constructs a single ``Vector`` from an existing ``Kokkos::View``. The View
      ExecutionSpace and MemorySpace must match the ExecutionSpace and
      MemorySpace provided as template arguments.

      :param view: A 1D ``Kokkos::View``
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: Vector(view_type view, host_view_type host_view, SUNContext sunctx)

      Constructs a single ``Vector`` from an existing ``Kokkos::View`` for the
      device and the host. The ExecutionSpace and MemorySpace of the device View
      must match the ExecutionSpace and MemorySpace provided as template arguments.

      :param view: A 1D ``Kokkos::View`` for the device
      :param host_view: A 1D ``Kokkos::View`` that is a ``Kokkos::HostMirrror`` for the device view
      :param sunctx: the SUNDIALS simulation context object (:c:type:`SUNContext`)

   .. cpp:function:: Vector(Vector&& that_vector) noexcept

      Move constructor.

   .. cpp:function:: Vector(const Vector& that_vector)

      Copy constructor. This creates a clone of the Vector, i.e., it creates
      a new Vector with the same properties, such as length, but it does not
      copy the data.

   .. cpp:function:: Vector& operator=(Vector&& rhs) noexcept

      Move assignment.

   .. cpp:function:: Vector& operator=(const Vector& rhs)

      Copy assignment. This creates a clone of the Vector, i.e., it creates
      a new Vector with the same properties, such as length, but it does not
      copy the data.

   .. cpp:function:: virtual ~Vector() = default;

      Default destructor.

   .. cpp:function:: size_type Length()

      Get the vector length i.e., ``extent(0)``.

   .. cpp:function:: view_type View()

      Get the underlying ``Kokkos:View`` for the device.

   .. cpp:function:: host_view_type HostView()

      Get the underlying ``Kokkos:View`` for the host.

   .. cpp:function:: operator N_Vector() override

      Implicit conversion to a :c:type:`N_Vector`.

   .. cpp:function:: operator N_Vector() const override

      Implicit conversion to a :c:type:`N_Vector`.

   .. cpp:function:: N_Vector get() override

      Explicit conversion to a :c:type:`N_Vector`.
      
      .. versionadded:: x.y.z
      
         Replaces the ``Convert`` method which was deprecated.

   .. cpp:function:: N_Vector get() const override

      Explicit conversion to a :c:type:`N_Vector`.
      
      .. versionadded:: x.y.z

         Replaces the ``Convert`` method which was deprecated.


.. cpp:function:: template<class VectorType> inline VectorType* GetVec(N_Vector v)

   Get the :cpp:type:`Vector` wrapped by a `N_Vector`.

.. cpp:function:: void CopyToDevice(N_Vector v)

   Copy the data from the host view to the device view with ``Kokkos::deep_copy``.

.. cpp:function:: void CopyFromDevice(N_Vector v)

   Copy the data to the host view from the device view with ``Kokkos::deep_copy``.

.. cpp:function:: template<class VectorType> void CopyToDevice(VectorType& v)

   Copy the data from the host view to the device view with ``Kokkos::deep_copy``.

.. cpp:function:: template<class VectorType> void CopyFromDevice(VectorType& v)

   Copy the data to the host view from the device view with ``Kokkos::deep_copy``.
