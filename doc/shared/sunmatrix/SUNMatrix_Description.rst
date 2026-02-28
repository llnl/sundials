..
   Programmer(s): Daniel R. Reynolds @ UMBC
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

.. _SUNMatrix.Description:

Description of the SUNMATRIX Modules
====================================

For problems that involve direct methods for solving linear systems,
the SUNDIALS packages not only operate on generic vectors, but also
on generic matrices (of type ``SUNMatrix``), through a set of
operations defined by the particular SUNMATRIX implementation.
Users can provide their own specific implementation of the
SUNMATRIX module, particularly in cases where they provide their
own ``N_Vector`` and/or linear solver modules, and require matrices
that are compatible with those implementations.  The generic
``SUNMatrix`` operations are described below, and descriptions of
the SUNMATRIX implementations provided with SUNDIALS follow.

The generic ``SUNMatrix`` type has been modeled after the
object-oriented style of the generic :c:type:`N_Vector` type.
Specifically, a generic ``SUNMatrix`` is a pointer to a structure
that has an implementation-dependent *content* field containing
the description and actual data of the matrix, and an *ops* field
pointing to a structure with generic matrix operations.

A :c:type:`SUNMatrix` is a pointer to the :c:struct:`_generic_SUNMatrix`
structure:

.. c:type:: struct _generic_SUNMatrix *SUNMatrix

.. c:struct:: _generic_SUNMatrix

   The structure defining the SUNDIALS matrix class.

   .. c:member:: void *content

      Pointer to matrix-specific member data

   .. c:member:: struct _generic_SUNMatrix_Ops *ops

      A virtual table of matrix operations provided by a specific
      implementation

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context

The virtual table structure is defined as

.. c:struct:: _generic_SUNMatrix_Ops

   The structure defining :c:type:`SUNMatrix` operations.

   .. c:member:: SUNMatrix_ID (*getid)(SUNMatrix)

      The function implementing :c:func:`SUNMatGetID`

   .. c:member:: SUNMatrix (*clone)(SUNMatrix)

      The function implementing :c:func:`SUNMatClone`

   .. c:member:: void (*destroy)(SUNMatrix)

      The function implementing :c:func:`SUNMatDestroy`

   .. c:member:: SUNErrCode (*zero)(SUNMatrix)

      The function implementing :c:func:`SUNMatZero`

   .. c:member:: SUNErrCode (*copy)(SUNMatrix, SUNMatrix)

      The function implementing :c:func:`SUNMatCopy`

   .. c:member:: SUNErrCode (*scaleadd)(sunrealtype, SUNMatrix, SUNMatrix)

      The function implementing :c:func:`SUNMatScaleAdd`

   .. c:member:: SUNErrCode (*scaleaddi)(sunrealtype, SUNMatrix)

      The function implementing :c:func:`SUNMatScaleAddI`

   .. c:member:: SUNErrCode (*matvecsetup)(SUNMatrix)

      The function implementing :c:func:`SUNMatMatvecSetup`

   .. c:member:: SUNErrCode (*matvec)(SUNMatrix, N_Vector, N_Vector)

      The function implementing :c:func:`SUNMatMatvec`

   .. c:member:: SUNErrCode (*mathermitiantransposevec)(SUNMatrix, N_Vector, N_Vector)

      The function implementing :c:func:`SUNMatHermitianTransposeVec`
      
      .. versionadded:: 7.3.0

   .. c:member:: SUNErrCode (*space)(SUNMatrix, long int*, long int*)

      The function implementing :c:func:`SUNMatSpace`


The generic SUNMATRIX module defines and implements the matrix
operations acting on a ``SUNMatrix``. These routines are nothing but
wrappers for the matrix operations defined by a particular SUNMATRIX
implementation, which are accessed through the *ops* field of the
``SUNMatrix`` structure. To illustrate this point we show below the
implementation of a typical matrix operation from the generic
SUNMATRIX module, namely ``SUNMatZero``, which sets all values of a
matrix ``A`` to zero, returning a flag denoting a successful/failed
operation:

.. code-block:: c

   SUNErrCode SUNMatZero(SUNMatrix A)
   {
     return(A->ops->zero(A));
   }

:numref:`SUNMatrix.Ops` contains a complete list of all
matrix operations defined by the generic SUNMATRIX module.  A
particular implementation of the SUNMATRIX module must:

* Specify the *content* field of the ``SUNMatrix`` object.

* Define and implement a minimal subset of the matrix operations.
  See the documentation for each SUNDIALS package and/or linear solver
  to determine which SUNMATRIX operations they require.

  Note that the names of these routines should be unique to that
  implementation in order to permit using more than one SUNMATRIX
  module (each with different ``SUNMatrix`` internal data
  representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free a ``SUNMatrix`` with the new *content*
  field and with *ops* pointing to the new matrix operations.

* Optionally, define and implement additional user-callable routines
  acting on the newly defined ``SUNMatrix`` (e.g., a routine to print the
  *content* for debugging purposes).

* Optionally, provide accessor macros as needed for that particular
  implementation to be used to access different parts in the content
  field of the newly defined ``SUNMatrix``.

To aid in the creation of custom SUNMATRIX modules the generic SUNMATRIX module
provides three utility functions :c:func:`SUNMatNewEmpty`,  :c:func:`SUNMatCopyOps()`,
and :c:func:`SUNMatFreeEmpty`. When used in custom SUNMATRIX constructors and clone
routines these functions will ease the introduction of any new optional matrix
operations to the SUNMATRIX API by ensuring only required operations need to be
set and all operations are copied when cloning a matrix.

.. c:function:: SUNMatrix SUNMatNewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNMatrix`` object and initializes its
  content pointer and the function pointers in the operations structure to ``NULL``.

  **Return value:**
     If successful, this function returns a ``SUNMatrix`` object. If an error
     occurs when allocating the object, then this routine will return ``NULL``.

.. c:function:: SUNErrCode SUNMatCopyOps(SUNMatrix A, SUNMatrix B)

  This function copies the function pointers in the ``ops`` structure of ``A``
  into the ``ops`` structure of ``B``.

   **Arguments:**
      * *A* -- the matrix to copy operations from.
      * *B* -- the matrix to copy operations to.

   **Return value:**
      * A :c:type:`SUNErrCode`

.. c:function:: void SUNMatFreeEmpty(SUNMatrix A)

  This routine frees the generic ``SUNMatrix`` object, under the assumption that any
  implementation-specific data that was allocated within the underlying content structure
  has already been freed. It will additionally test whether the ops pointer is ``NULL``,
  and, if it is not, it will free it as well.

   **Arguments:**
      * *A* -- the SUNMatrix object to free


.. c:type:: SUNMatrix_ID

   Each SUNMATRIX implementation included in SUNDIALS has a unique identifier
   specified in enumeration and shown in
   :numref:`SUNMatrix.Description.matrixIDs`. It is recommended that a
   user-supplied SUNMATRIX implementation use the ``SUNMATRIX_CUSTOM``
   identifier.


.. _SUNMatrix.Description.matrixIDs:
.. table:: Identifiers associated with matrix kernels supplied with SUNDIALS
   :align: center

   ======================  =================================================
   Matrix ID               Matrix type
   ======================  =================================================
   SUNMATRIX_BAND          Band :math:`M \times M` matrix
   SUNMATRIX_CUSPARSE      CUDA sparse CSR matrix
   SUNMATRIX_CUSTOM        User-provided custom matrix
   SUNMATRIX_DENSE         Dense :math:`M \times N` matrix
   SUNMATRIX_GINKGO        SUNMatrix wrapper for Ginkgo matrices
   SUNMATRIX_MAGMADENSE    Dense :math:`M \times N` matrix
   SUNMATRIX_ONEMKLDENSE   oneMKL dense :math:`M \times N` matrix
   SUNMATRIX_SLUNRLOC      SUNMatrix wrapper for SuperLU_DIST SuperMatrix
   SUNMATRIX_SPARSE        Sparse (CSR or CSC) :math:`M\times N` matrix
   ======================  =================================================
