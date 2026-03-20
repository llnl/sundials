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

.. _SUNMemory.CUDA:

The SUNMemoryHelper_Cuda Implementation
=======================================

The SUNMemoryHelper_Cuda module is an implementation of the ``SUNMemoryHelper``
API that interfaces to the NVIDIA :cite:p:`cuda_site` library.  The
implementation defines the constructor

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Cuda(SUNContext sunctx)

   Allocates and returns a ``SUNMemoryHelper`` object for handling CUDA memory
   if successful. Otherwise it returns ``NULL``.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNMemoryHelper` object; otherwise it will return ``NULL``.


.. _SUNMemory.CUDA.Operations:

SUNMemoryHelper_Cuda API Functions
----------------------------------

The implementation provides the following operations defined by the
``SUNMemoryHelper`` API:

.. c:function:: SUNMemory SUNMemoryHelper_Alloc_Cuda(SUNMemoryHelper helper, \
                                                     SUNMemory memptr, \
                                                     size_t mem_size, \
                                                     SUNMemoryType mem_type, \
                                                     void* queue)

   Allocates a ``SUNMemory`` object whose ``ptr`` field is allocated for
   ``mem_size`` bytes and is of type ``mem_type``. The new object will have
   ownership of ``ptr`` and will be deallocated when
   :c:func:`SUNMemoryHelper_Dealloc` is called.

   :param helper: the ``SUNMemoryHelper`` object.
   :param memptr: pointer to the allocated ``SUNMemory``.
   :param mem_size: the size in bytes of the ``ptr``.
   :param mem_type: the ``SUNMemoryType`` of the ``ptr``. Supported values are:
     * ``SUNMEMTYPE_HOST`` -- memory is allocated with a call to ``malloc``.
     * ``SUNMEMTYPE_PINNED`` -- memory is allocated with a call to ``cudaMallocHost``.
     * ``SUNMEMTYPE_DEVICE`` -- memory is allocated with a call to ``cudaMalloc``.
     * ``SUNMEMTYPE_UVM`` -- memory is allocated with a call to ``cudaMallocManaged``.
   :param queue: currently unused.
   :return: A new :c:type:`SUNMemory` object.


.. c:function:: SUNMemory SUNMemoryHelper_AllocStrided_Cuda(SUNMemoryHelper helper, \
                                                           SUNMemory memptr, \
                                                           size_t mem_size, \
                                                           size_t stride, \
                                                           SUNMemoryType mem_type, \
                                                           void* queue)

    Allocates a ``SUNMemory`` object whose ``ptr`` field is allocated for
    ``mem_size`` bytes and is of type ``mem_type``. The new object will have 
    ownership of ``ptr`` and will be deallocated when :c:func:`SUNMemoryHelper_Dealloc` is called.

    :param helper: the ``SUNMemoryHelper`` object.
    :param memptr: pointer to the allocated ``SUNMemory``.
    :param mem_size: the size in bytes of the ``ptr``.
    :param stride: the number of bytes between elements in the array.
    :param mem_type: the ``SUNMemoryType`` of the ``ptr``.
    :param queue: currently unused.
    :return: A new :c:type:`SUNMemory` object.
    
    .. versionadded:: 7.3.0
    

.. c:function:: SUNErrCode SUNMemoryHelper_Dealloc_Cuda(SUNMemoryHelper helper, \
                                                 SUNMemory mem, void* queue)

   Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
   deallocates the ``mem`` object.

   :param helper: the ``SUNMemoryHelper`` object.
   :param mem: the ``SUNMemory`` object.
   :param queue: currently unused.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNMemoryHelper_Copy_Cuda(SUNMemoryHelper helper, \
                                              SUNMemory dst, SUNMemory src, \
                                              size_t mem_size, void* queue)

   Synchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper`` object
   will use the memory types of ``dst`` and ``src`` to determine the appropriate
   transfer type necessary.

   :param helper: the ``SUNMemoryHelper`` object.
   :param dst: the destination memory to copy to.
   :param src: the source memory to copy from.
   :param mem_size: the number of bytes to copy.
   :param queue: currently unused.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNMemoryHelper_CopyAsync_Cuda(SUNMemoryHelper helper, \
                                                   SUNMemory dst, \
                                                   SUNMemory src, \
                                                   size_t mem_size, void* queue)

   Asynchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper`` object
   will use the memory types of ``dst`` and ``src`` to determine the appropriate
   transfer type necessary.

   :param helper: the ``SUNMemoryHelper`` object.
   :param dst: the destination memory to copy to.
   :param src: the source memory to copy from.
   :param mem_size: the number of bytes to copy.
   :param queue: the ``cudaStream_t`` handle for the stream that the copy will
     be performed on.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNMemoryHelper_GetAllocStats_Cuda(SUNMemoryHelper helper, SUNMemoryType mem_type, unsigned long* num_allocations, \
                                                       unsigned long* num_deallocations, size_t* bytes_allocated, \
                                                       size_t* bytes_high_watermark)

   Returns statistics about memory allocations performed with the helper.

   :param helper: the ``SUNMemoryHelper`` object.
   :param mem_type: the ``SUNMemoryType`` to get stats for.
   :param num_allocations: (output argument) number of memory allocations done through the helper.
   :param num_deallocations: (output argument) number of memory deallocations done through the helper.
   :param bytes_allocated: (output argument) total number of bytes allocated through the helper at the moment this function is called.
   :param bytes_high_watermark: (output argument) max number of bytes allocated through the helper at any moment in the lifetime of the helper.
   :return: A :c:type:`SUNErrCode` indicating success or failure.
