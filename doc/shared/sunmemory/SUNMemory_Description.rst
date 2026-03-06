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

.. _SUNMemory.Description:

The SUNMemoryHelper API
=======================

This API consists of three new SUNDIALS types: :c:type:`SUNMemoryType`,
:c:type:`SUNMemory`, and :c:type:`SUNMemoryHelper`:


.. c:type:: struct SUNMemory_ *SUNMemory

   The :c:type:`SUNMemory` type is a pointer the structure

   .. c:struct:: SUNMemory_

      .. c:member:: void* ptr;

         The  actual data.

      .. c:member:: SUNMemoryType type;

         The data memory type.

      .. c:member:: sunbooleantype own;

         A flag indicating ownership.

      .. c:member:: size_t bytes;

         The size of the data allocated.

      .. c:member:: size_t stride;

         .. versionadded:: 7.3.0

         The stride of the data.

.. c:function:: SUNMemory SUNMemoryNewEmpty(SUNContext sunctx)

   This function returns an empty :c:type:`SUNMemory` object.

   :param sunctx: the :c:type:`SUNContext` object.

   :return: an uninitialized :c:type:`SUNMemory` object

   .. versionchanged:: 7.0.0

      The function signature was updated to add the :c:type:`SUNContext` argument.


.. c:enum:: SUNMemoryType

   The :c:type:`SUNMemoryType` type is an enumeration that defines the supported
   memory types:

   .. c:enumerator:: SUNMEMTYPE_HOST

      Pageable memory accessible on the host

   .. c:enumerator:: SUNMEMTYPE_PINNED

      Page-locked memory accessible on the host

   .. c:enumerator:: SUNMEMTYPE_DEVICE

      Memory accessible from the device

   .. c:enumerator:: SUNMEMTYPE_UVM

      Memory accessible from the host or device


.. c:type:: struct SUNMemoryHelper_ *SUNMemoryHelper

   The :c:type:`SUNMemoryHelper` type is a pointer to the structure

   .. c:struct:: SUNMemoryHelper_

      .. c:member:: void* content;

         Pointer to the implementation-specific member data

      .. c:member:: void* queue;

         Pointer to the implementation-specific queue (e.g., a ``cudaStream_t*``) 
         to use by default when one is not provided for an operation
         
         .. versionadded:: 7.3.0

      .. c:member:: SUNMemoryHelper_Ops ops;

         A virtual method table of member functions

      .. c:member:: SUNContext sunctx;

         The SUNDIALS simulation context


.. c:type:: struct SUNMemoryHelper_Ops_ *SUNMemoryHelper_Ops

   The ``SUNMemoryHelper_Ops`` type is defined as a pointer to the structure
   containing the function pointers to the member function implementations. This
   structure is define as

   .. c:struct:: SUNMemoryHelper_Ops_

      .. c:member:: SUNErrCode (*alloc)(SUNMemoryHelper, SUNMemory* memptr, size_t mem_size, SUNMemoryType mem_type, void* queue)

         The function implementing :c:func:`SUNMemoryHelper_Alloc`

      .. c:member:: SUNErrCode (*allocstrided)(SUNMemoryHelper, SUNMemory* memptr, size_t mem_size, size_t stride, SUNMemoryType mem_type, void* queue)

         The function implementing :c:func:`SUNMemoryHelper_AllocStrided`
         
         .. versionadded:: 7.3.0

      .. c:member:: SUNErrCode (*dealloc)(SUNMemoryHelper, SUNMemory mem, void* queue)

         The function implementing :c:func:`SUNMemoryHelper_Dealloc`

      .. c:member:: SUNErrCode (*copy)(SUNMemoryHelper, SUNMemory dst, SUNMemory src, size_t mem_size, void* queue)

         The function implementing :c:func:`SUNMemoryHelper_Copy`

      .. c:member:: SUNErrCode (*copyasync)(SUNMemoryHelper, SUNMemory dst, SUNMemory src, size_t mem_size, void* queue)

         The function implementing :c:func:`SUNMemoryHelper_CopyAsync`

      .. c:member:: SUNErrCode (*getallocstats)(SUNMemoryHelper, SUNMemoryType mem_type, unsigned long* num_allocations, unsigned long* num_deallocations, size_t* bytes_allocated, size_t* bytes_high_watermark)

         The function implementing :c:func:`SUNMemoryHelper_GetAllocStats`

      .. c:member:: SUNMemoryHelper (*clone)(SUNMemoryHelper)

         The function implementing :c:func:`SUNMemoryHelper_Clone`

      .. c:member:: SUNErrCode (*destroy)(SUNMemoryHelper)

         The function implementing :c:func:`SUNMemoryHelper_Destroy`


.. _SUNMemory.Description.Required:

Implementation defined operations
---------------------------------

The SUNMemory API defines the following operations that an implementation to
must define:

.. c:function:: SUNMemory SUNMemoryHelper_Alloc(SUNMemoryHelper helper, \
                                                SUNMemory* memptr, \
                                                size_t mem_size, \
                                                SUNMemoryType mem_type, \
                                                void* queue)

   Allocates a :c:type:`SUNMemory` object whose ``ptr`` field is allocated for
   ``mem_size`` bytes and is of type ``mem_type``. The new object will have
   ownership of ``ptr`` and will be deallocated when
   :c:func:`SUNMemoryHelper_Dealloc` is called.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param memptr: pointer to the allocated :c:type:`SUNMemory`.
   :param mem_size: the size in bytes of the ``ptr``.
   :param mem_type: the :c:type:`SUNMemoryType` of the ``ptr``.
   :param queue: typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   :return: A new :c:type:`SUNMemory` object


.. c:function:: SUNMemory SUNMemoryHelper_AllocStrided(SUNMemoryHelper helper, \
                                                       SUNMemory* memptr, \
                                                       size_t mem_size, size_t stride, \
                                                       SUNMemoryType mem_type, \
                                                       void* queue)

   Allocates a :c:type:`SUNMemory` object whose ``ptr`` field is allocated for
   ``mem_size`` bytes with the specified stride, and is of type ``mem_type``.
   The new object will have ownership of ``ptr`` and will be deallocated when
   :c:func:`SUNMemoryHelper_Dealloc` is called.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param memptr: pointer to the allocated :c:type:`SUNMemory`.
   :param mem_size: the size in bytes of the ``ptr``.
   :param stride: the stride of the memory in bytes.
   :param mem_type: the :c:type:`SUNMemoryType` of the ``ptr``.
   :param queue: typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   :return: A new :c:type:`SUNMemory` object
   
   .. versionadded:: 7.3.0


.. c:function:: SUNErrCode SUNMemoryHelper_Dealloc(SUNMemoryHelper helper, \
                                            SUNMemory mem, void* queue)

   Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
   deallocates the ``mem`` object.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param mem: the :c:type:`SUNMemory` object.
   :param queue: typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNMemoryHelper_Copy(SUNMemoryHelper helper, \
                                         SUNMemory dst, SUNMemory src, \
                                         size_t mem_size, void* queue)

   Synchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper``
   object should use the memory types of ``dst`` and ``src`` to determine
   the appropriate transfer type necessary.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param dst: the destination memory to copy to.
   :param src: the source memory to copy from.
   :param mem_size: the number of bytes to copy.
   :param queue: typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   :return: A :c:type:`SUNErrCode` indicating success or failure.



.. _SUNMemory.Description.Utilities:

Utility Functions
-----------------

The SUNMemoryHelper API defines the following functions which do not
require a SUNMemoryHelper instance:

.. c:function:: SUNMemory SUNMemoryHelper_Alias(SUNMemoryHelper helper, SUNMemory mem1)

   Returns a :c:type:`SUNMemory` object whose ``ptr`` field points to the same address
   as ``mem1``. The new object *will not* have ownership of ``ptr``, therefore,
   it will not free ``ptr`` when :c:func:`SUNMemoryHelper_Dealloc` is called.

   :param helper: a :c:type:`SUNMemoryHelper` object.
   :param mem1: a :c:type:`SUNMemory` object.


   :return: A :c:type:`SUNMemory` object or ``NULL`` if an error occurs.

   .. versionchanged:: 7.0.0

      The :c:type:`SUNMemoryHelper` argument was added to the function signature.


.. c:function:: SUNMemory SUNMemoryHelper_Wrap(SUNMemoryHelper helper, void* ptr, \
                                               SUNMemoryType mem_type)

   Returns a :c:type:`SUNMemory` object whose ``ptr`` field points to the ``ptr``
   argument passed to the function. The new object *will not* have ownership of
   ``ptr``, therefore, it will not free ``ptr`` when
   :c:func:`SUNMemoryHelper_Dealloc` is called.

   :param helper: a :c:type:`SUNMemoryHelper` object.
   :param ptr: the data pointer to wrap in a :c:type:`SUNMemory` object.
   :param mem_type: the :c:type:`SUNMemoryType` of the ``ptr``.


   :return: A :c:type:`SUNMemory` object or ``NULL`` if an error occurs.

   .. versionchanged:: 7.0.0

      The :c:type:`SUNMemoryHelper` argument was added to the function signature.


.. c:function:: SUNMemoryHelper SUNMemoryHelper_NewEmpty(SUNContext sunctx)

   Returns an empty :c:type:`SUNMemoryHelper`. This is useful for building custom
   :c:type:`SUNMemoryHelper` implementations.

   :param helper: a :c:type:`SUNMemoryHelper` object.


   :return: A :c:type:`SUNMemoryHelper` object or ``NULL`` if an error occurs.

   .. versionchanged:: 7.0.0

      The :c:type:`SUNMemoryHelper` argument was added to the function signature.


.. c:function:: SUNErrCode SUNMemoryHelper_CopyOps(SUNMemoryHelper src, \
                                            SUNMemoryHelper dst)

   Copies the ``ops`` field of ``src`` to the ``ops`` field of ``dst``.
   This is useful for building custom :c:type:`SUNMemoryHelper` implementations.

   :param src: the object to copy from.
   :param dst: the object to copy to.


   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNMemoryHelper_GetAllocStats(SUNMemoryHelper helper, SUNMemoryType mem_type, unsigned long* num_allocations, \
                                                  unsigned long* num_deallocations, size_t* bytes_allocated, \
                                                  size_t* bytes_high_watermark)

   Returns statistics about the allocations performed with the helper.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param mem_type: the :c:type:`SUNMemoryType` to get stats for.
   :param num_allocations:  (output argument) number of allocations done through the helper.
   :param num_deallocations:  (output argument) number of deallocations done through the helper.
   :param bytes_allocated:  (output argument) total number of bytes allocated through the helper at the moment this function is called.
   :param bytes_high_watermark:  (output argument) max number of bytes allocated through the helper at any moment in the lifetime of the helper.


   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNMemoryHelper_SetDefaultQueue(SUNMemoryHelper helper, void* queue)

   Sets the default queue for the helper.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param queue: pointer to the queue to use by default.
   :return: A :c:type:`SUNErrCode` indicating success or failure.
   
   .. versionadded:: 7.3.0


.. _SUNMemory.Description.Overridable:

Implementation overridable operations with defaults
---------------------------------------------------

In addition, the SUNMemoryHelper API defines the following *optionally
overridable* operations which an implementation may define:


.. c:function:: SUNErrCode SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, \
                                              SUNMemory dst, SUNMemory src, \
                                              size_t mem_size, void* queue)

   Asynchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper`` object
   should use the memory types of ``dst`` and ``src`` to determine the
   appropriate transfer type necessary.  The ``ctx`` argument is used when a
   different execution stream needs to be provided to perform the copy in,
   e.g. with ``CUDA`` this would be a ``cudaStream_t``.

   :param helper: the :c:type:`SUNMemoryHelper` object.
   :param dst: the destination memory to copy to.
   :param src: the source memory to copy from.
   :param mem_size: the number of bytes to copy.
   :param queue: typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.


   An ``int`` flag indicating success (zero) or failure (non-zero).

   .. note::

      If this operation is not defined by the implementation, then
      :c:func:`SUNMemoryHelper_Copy` will be used.

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper helper)

   Clones the :c:type:`SUNMemoryHelper` object itself.

   :param helper: the :c:type:`SUNMemoryHelper` object to clone.

   :return: A :c:type:`SUNMemoryHelper` object.

   .. note::

      If this operation is not defined by the implementation, then the default
      clone will only copy the ``SUNMemoryHelper_Ops`` structure stored in
      ``helper->ops``, and not the ``helper->content`` field.


.. c:function:: SUNErrCode SUNMemoryHelper_Destroy(SUNMemoryHelper helper)

   Destroys (frees) the :c:type:`SUNMemoryHelper` object itself.

   :param helper: the :c:type:`SUNMemoryHelper` object to destroy.

   :return: A :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      If this operation is not defined by the implementation, then the default
      destroy will only free the ``helper->ops`` field and the ``helper``
      itself. The ``helper->content`` field will not be freed.


.. _SUNMemory.Description.Custom:

Implementing a custom SUNMemoryHelper
-------------------------------------

A particular implementation of the SUNMemoryHelper API must:

*  Define and implement the required operations. Note that the names of
   these routines should be unique to that implementation in order to
   permit using more than one SUNMemoryHelper module in the same code.

*  Optionally, specify the *content* field of SUNMemoryHelper.

*  Optionally, define and implement additional user-callable routines
   acting on the newly defined SUNMemoryHelper.

An example of a custom SUNMemoryHelper is given in
``examples/utilities/custom_memory_helper.h``.
