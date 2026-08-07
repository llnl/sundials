.. ----------------------------------------------------------------
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

.. _SUNAdjoint.CheckpointScheme:

The SUNAdjointCheckpointScheme Class
====================================

.. versionadded:: 7.3.0

As with other SUNDIALS classes, the :c:type:`SUNAdjointCheckpointScheme` abstract base class is
implemented using a C structure containing a ``content`` pointer to the derived class member data
and a structure of function pointers to the derived class implementations of the virtual methods.

.. c:type:: SUNAdjointCheckpointScheme

   A class that provides an interface for checkpointing states during forward
   integration and accessing them as needed during the backwards integration
   of the adjoint model.

.. c:enum:: SUNDataIOMode

   .. c:enumerator:: SUNDATAIOMODE_INMEM

      The IO mode for data that is stored in addressable random access memory.
      The location of the memory (e.g., CPU or GPU) is not specified by this mode.


.. _SUNAdjoint.CheckpointScheme.BaseClassMethods:

Base Class Methods
------------------

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx, \
   SUNAdjointCheckpointScheme* cs_ptr)

   :param sunctx: The SUNDIALS simulation context
   :param cs_ptr: on output, a pointer to a new :c:type:`SUNAdjointCheckpointScheme` object

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving(SUNAdjointCheckpointScheme self, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   Determines if the (step_num, stage_num) should be checkpointed or not.

   :param self: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param yes_or_no: boolean indicating if the checkpoint should be saved or not

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme self, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, N_Vector y)

   Inserts the vector as the checkpoint for (step_num, stage_num).

   :param self: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the time of the checkpoint
   :param y: the state vector to checkpoint

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_LoadVector(SUNAdjointCheckpointScheme self, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype peek, N_Vector* yout, sunrealtype* tout)

   Loads the checkpointed vector for (step_num, stage_num).

   :param self: the :c:type:`SUNAdjointCheckpointScheme` object
   :param step_num: the step number of the checkpoint
   :param stage_num: the stage number of the checkpoint
   :param t: the desired time of the checkpoint
   :param peek: if true, then the checkpoint will be loaded but not deleted regardless
      of other implementation-specific settings. If false, then the checkpoint may be
      deleted depending on the implementation.
   :param yout: the loaded state vector
   :param tout: on output, the time of the checkpoint

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme self, \
   sunbooleantype on_or_off)

   Enables or disables dense checkpointing (checkpointing every step/stage). When dense checkpointing
   is disabled, the checkpointing interval that was set when the object was created is restored.

   :param self: the :c:type:`SUNAdjointCheckpointScheme` object
   :param on_or_off: if true, dense checkpointing will be turned on, if false it will be turned off.

   :returns: A :c:type:`SUNErrCode` indicating failure or success.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme* cs_ptr)

   Destroys (deallocates) the SUNAdjointCheckpointScheme object.

   :param cs_ptr: pointer to a :c:type:`SUNAdjointCheckpointScheme` object

   :returns: A :c:type:`SUNErrCode` indicating failure or success.


.. _SUNAdjoint.CheckpointScheme.ImplMethods:

Implementation Specific Methods
-------------------------------

This section describes the virtual methods defined by the :c:type:`SUNAdjointCheckpointScheme`
abstract base class.

.. c:type:: SUNErrCode (*SUNAdjointCheckpointSchemeNeedsSavingFn)(SUNAdjointCheckpointScheme check_scheme, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no)

   This type represents a function with the signature of
   :c:func:`SUNAdjointCheckpointScheme_NeedsSaving`.

.. c:type:: SUNErrCode (*SUNAdjointCheckpointSchemeInsertVectorFn)(SUNAdjointCheckpointScheme check_scheme, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, N_Vector y)

   This type represents a function with the signature of
   :c:func:`SUNAdjointCheckpointScheme_InsertVector`.

.. c:type:: SUNErrCode (*SUNAdjointCheckpointSchemeLoadVectorFn)(SUNAdjointCheckpointScheme check_scheme, \
   suncountertype step_num, suncountertype stage_num, sunrealtype t, sunbooleantype peek, N_Vector* yout, sunrealtype* tout)

   This type represents a function with the signature of
   :c:func:`SUNAdjointCheckpointScheme_LoadVector`.

.. c:type:: SUNErrCode (*SUNAdjointCheckpointSchemeEnableDenseFn)(SUNAdjointCheckpointScheme check_scheme, \
   sunbooleantype on_or_off)

   This type represents a function with the signature of
   :c:func:`SUNAdjointCheckpointScheme_EnableDense`.

.. c:type:: SUNErrCode (*SUNAdjointCheckpointSchemeDestroyFn)(SUNAdjointCheckpointScheme* check_scheme_ptr)

   This type represents a function with the signature of
   :c:func:`SUNAdjointCheckpointScheme_Destroy`.


.. _SUNAdjoint.CheckpointScheme.SetContentMembers:

Setting Content and Member Functions
------------------------------------

These functions can be used to set the content pointer or virtual method pointers
as needed when implementing the abstract base class.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_SetNeedsSavingFn(SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeNeedsSavingFn fn)

   This function attaches a :c:type:`SUNAdjointCheckpointSchemeNeedsSavingFn` function to a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param fn: the :c:type:`SUNAdjointCheckpointSchemeNeedsSavingFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_SetInsertVectorFn(SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeInsertVectorFn fn)

   This function attaches a :c:type:`SUNAdjointCheckpointSchemeInsertVectorFn` function to a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param fn: the :c:type:`SUNAdjointCheckpointSchemeInsertVectorFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_SetLoadVectorFn(SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeLoadVectorFn fn)

   This function attaches a :c:type:`SUNAdjointCheckpointSchemeLoadVectorFn` function to a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param fn: the :c:type:`SUNAdjointCheckpointSchemeLoadVectorFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_SetDestroyFn(SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeDestroyFn fn)

   This function attaches a :c:type:`SUNAdjointCheckpointSchemeDestroyFn` function to a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param fn: the :c:type:`SUNAdjointCheckpointSchemeDestroyFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_SetEnableDenseFn(SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeEnableDenseFn fn)

   This function attaches a :c:type:`SUNAdjointCheckpointSchemeEnableDenseFn` function to a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param fn: the :c:type:`SUNAdjointCheckpointSchemeEnableDenseFn` function to attach.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_SetContent(SUNAdjointCheckpointScheme self, void* content)

   This function attaches a member data (content) pointer to a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param content: a pointer to the checkpoint scheme member data.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_GetContent(SUNAdjointCheckpointScheme self, void** content)

   This function retrieves the member data (content) pointer from a
   :c:type:`SUNAdjointCheckpointScheme` object.

   :param self: a checkpoint scheme object.
   :param content: a pointer to set to the checkpoint scheme member data pointer.
   :return: A :c:type:`SUNErrCode` indicating success or failure.


.. _SUNAdjoint.CheckpointScheme.Fixed:

The SUNAdjointCheckpointScheme_Fixed Module
===========================================

The ``SUNAdjointCheckpointScheme_Fixed`` module implements a scheme where a checkpoint is saved at some
fixed interval (in time steps). The module supports checkpointing of time step states only, or time step
states with intermediate stage states as well (for multistage methods). When used with a
fixed time step size then the number of checkpoints that will be saved is fixed. However, with
adaptive time steps the number of checkpoints stored with this scheme is unbounded.

The diagram below illustrates how checkpoints are stored with this scheme:

.. figure:: /figs/sunadjoint_ckpt_fixed.png
   :width: 75 %
   :align: center


Base-class Method Overrides
---------------------------

The ``SUNAdjointCheckpointScheme_Fixed`` module implements the following :c:type:`SUNAdjointCheckpointScheme` functions:

* :c:func:`SUNAdjointCheckpointScheme_NeedsSaving`
* :c:func:`SUNAdjointCheckpointScheme_InsertVector`
* :c:func:`SUNAdjointCheckpointScheme_LoadVector`
* :c:func:`SUNAdjointCheckpointScheme_Destroy`
* :c:func:`SUNAdjointCheckpointScheme_EnableDense`


Implementation Specific Methods
-------------------------------

The ``SUNAdjointCheckpointScheme_Fixed`` module also implements the following module-specific functions:

.. c:function:: SUNErrCode SUNAdjointCheckpointScheme_Create_Fixed(SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval, suncountertype estimate, sunbooleantype keep, SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)

   Creates a new :c:type:`SUNAdjointCheckpointScheme` object that checkpoints at a fixed interval.

   :param io_mode: The IO mode used for storing the checkpoints.
   :param mem_helper: Memory helper for managing memory.
   :param interval: The interval (in steps) between checkpoints.
   :param estimate: An estimate of the total number of checkpoints needed.
   :param keep: Keep data stored even after it is not needed anymore.
   :param sunctx: The :c:type:`SUNContext` for the simulation.
   :param check_scheme_ptr: Pointer to the newly constructed object.
   :returns: A :c:type:`SUNErrCode` indicating success or failure.
