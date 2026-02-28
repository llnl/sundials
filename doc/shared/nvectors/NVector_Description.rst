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

.. _NVectors.Description:

Description of the NVECTOR Modules
==================================

SUNDIALS solvers are written in a data-independent manner. They all operate on
generic vectors (of type :c:type:`N_Vector`) through a set of operations defined
by, and specific to, the particular vector implementation. Users can provide a
custom vector implementation or use one provided with SUNDIALS. The generic
operations are described below. In the sections following, the implementations
provided with SUNDIALS are described.

An :c:type:`N_Vector` is a pointer to the :c:struct:`_generic_N_Vector`
structure:

.. c:type:: struct _generic_N_Vector *N_Vector

.. c:struct:: _generic_N_Vector

   The structure defining the SUNDIALS vector class.

   .. c:member:: void *content

      Pointer to vector-specific member data.

   .. c:member:: N_Vector_Ops ops

      A virtual table of vector operations provided by a specific
      implementation.

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context

The virtual table structure is defined as

.. c:type:: _generic_N_Vector_Ops *N_Vector_Ops

.. c:struct:: _generic_N_Vector_Ops

   The structure defining :c:type:`N_Vector` operations.

   .. c:member:: N_Vector_ID (*nvgetvectorid)(N_Vector)

      The function implementing :c:func:`N_VGetVectorID`

   .. c:member:: N_Vector (*nvclone)(N_Vector)

      The function implementing :c:func:`N_VClone`

   .. c:member:: N_Vector (*nvcloneempty)(N_Vector)

      The function implementing :c:func:`N_VCloneEmpty`

   .. c:member:: void (*nvdestroy)(N_Vector)

      The function implementing :c:func:`N_VDestroy`

   .. c:member:: void (*nvspace)(N_Vector, sunindextype*, sunindextype*)

      The function implementing :c:func:`N_VSpace`

   .. c:member:: sunrealtype* (*nvgetarraypointer)(N_Vector)

      The function implementing :c:func:`N_VGetArrayPointer`

   .. c:member:: sunrealtype* (*nvgetdevicearraypointer)(N_Vector)

      The function implementing :c:func:`N_VGetDeviceArrayPointer`

   .. c:member:: void (*nvsetarraypointer)(sunrealtype*, N_Vector)

      The function implementing :c:func:`N_VSetArrayPointer`

   .. c:member:: SUNComm (*nvgetcommunicator)(N_Vector)

      The function implementing :c:func:`N_VGetCommunicator`

   .. c:member:: sunindextype (*nvgetlength)(N_Vector)

      The function implementing :c:func:`N_VGetLength`

   .. c:member:: sunindextype (*nvgetlocallength)(N_Vector)

      The function implementing :c:func:`N_VGetLocalLength`

   .. c:member:: void (*nvlinearsum)(sunrealtype, N_Vector, sunrealtype, N_Vector, N_Vector)

      The function implementing :c:func:`N_VLinearSum`

   .. c:member:: void (*nvconst)(sunrealtype, N_Vector)

      The function implementing :c:func:`N_VConst`

   .. c:member:: void (*nvprod)(N_Vector, N_Vector, N_Vector)

      The function implementing :c:func:`N_VProd`

   .. c:member:: void (*nvdiv)(N_Vector, N_Vector, N_Vector)

      The function implementing :c:func:`N_VDiv`

   .. c:member:: void (*nvscale)(sunrealtype, N_Vector, N_Vector)

      The function implementing :c:func:`N_VScale`

   .. c:member:: void (*nvabs)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VAbs`

   .. c:member:: void (*nvinv)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VInv`

   .. c:member:: void (*nvaddconst)(N_Vector, sunrealtype, N_Vector)

      The function implementing :c:func:`N_VAddConst`

   .. c:member:: sunrealtype (*nvdotprod)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VDotProd`

   .. c:member:: sunrealtype (*nvmaxnorm)(N_Vector)

      The function implementing :c:func:`N_VMaxNorm`

   .. c:member:: sunrealtype (*nvwrmsnorm)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VWrmsNorm`

   .. c:member:: sunrealtype (*nvwrmsnormmask)(N_Vector, N_Vector, N_Vector)

      The function implementing :c:func:`N_VWrmsNormMask`

   .. c:member:: sunrealtype (*nvmin)(N_Vector)

      The function implementing :c:func:`N_VMin`

   .. c:member:: sunrealtype (*nvwl2norm)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VWL2Norm`

   .. c:member:: sunrealtype (*nvl1norm)(N_Vector)

      The function implementing :c:func:`N_VL1Norm`

   .. c:member:: void (*nvcompare)(sunrealtype, N_Vector, N_Vector)

      The function implementing :c:func:`N_VCompare`

   .. c:member:: sunbooleantype (*nvinvtest)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VInvTest`

   .. c:member:: sunbooleantype (*nvconstrmask)(N_Vector, N_Vector, N_Vector)

      The function implementing :c:func:`N_VConstrMask`

   .. c:member:: sunrealtype (*nvminquotient)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VMinQuotient`

   .. c:member:: SUNErrCode (*nvlinearcombination)(int, sunrealtype*, N_Vector*, N_Vector)

      The function implementing :c:func:`N_VLinearCombination`

   .. c:member:: SUNErrCode (*nvscaleaddmulti)(int, sunrealtype*, N_Vector, N_Vector*, N_Vector*)

      The function implementing :c:func:`N_VScaleAddMulti`

   .. c:member:: SUNErrCode (*nvdotprodmulti)(int, N_Vector, N_Vector*, sunrealtype*)

      The function implementing :c:func:`N_VDotProdMulti`

   .. c:member:: SUNErrCode (*nvlinearsumvectorarray)(int, sunrealtype, N_Vector*, sunrealtype, N_Vector*, N_Vector*)

      The function implementing :c:func:`N_VLinearSumVectorArray`

   .. c:member:: SUNErrCode (*nvscalevectorarray)(int, sunrealtype*, N_Vector*, N_Vector*)

      The function implementing :c:func:`N_VScaleVectorArray`

   .. c:member:: SUNErrCode (*nvconstvectorarray)(int, sunrealtype, N_Vector*)

      The function implementing :c:func:`N_VConstVectorArray`

   .. c:member:: SUNErrCode (*nvwrmsnormvectorarray)(int, N_Vector*, N_Vector*, sunrealtype*)

      The function implementing :c:func:`N_VWrmsNormVectorArray`

   .. c:member:: SUNErrCode (*nvwrmsnormmaskvectorarray)(int, N_Vector*, N_Vector*, N_Vector, sunrealtype*)

      The function implementing :c:func:`N_VWrmsNormMaskVectorArray`

   .. c:member:: SUNErrCode (*nvscaleaddmultivectorarray)(int, int, sunrealtype*, N_Vector*, N_Vector**, N_Vector**)

      The function implementing :c:func:`N_VScaleAddMultiVectorArray`

   .. c:member:: SUNErrCode (*nvlinearcombinationvectorarray)(int, int, sunrealtype*, N_Vector**, N_Vector*)

      The function implementing :c:func:`N_VLinearCombinationVectorArray`

   .. c:member:: sunrealtype (*nvdotprodlocal)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VDotProdLocal`

   .. c:member:: sunrealtype (*nvmaxnormlocal)(N_Vector)

      The function implementing :c:func:`N_VMaxNormLocal`

   .. c:member:: sunrealtype (*nvminlocal)(N_Vector)

      The function implementing :c:func:`N_VMinLocal`

   .. c:member:: sunrealtype (*nvl1normlocal)(N_Vector)

      The function implementing :c:func:`N_VL1NormLocal`

   .. c:member:: sunbooleantype (*nvinvtestlocal)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VInvTestLocal`

   .. c:member:: sunbooleantype (*nvconstrmasklocal)(N_Vector, N_Vector, N_Vector)

      The function implementing :c:func:`N_VConstrMaskLocal`

   .. c:member:: sunrealtype (*nvminquotientlocal)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VMinQuotientLocal`

   .. c:member:: sunrealtype (*nvwsqrsumlocal)(N_Vector, N_Vector)

      The function implementing :c:func:`N_VWSqrSumLocal`

   .. c:member:: sunrealtype (*nvwsqrsummasklocal)(N_Vector, N_Vector, N_Vector)

      The function implementing :c:func:`N_VWSqrSumMaskLocal`

   .. c:member:: SUNErrCode (*nvdotprodmultilocal)(int, N_Vector, N_Vector*, sunrealtype*)

      The function implementing :c:func:`N_VDotProdMultiLocal`

   .. c:member:: SUNErrCode (*nvdotprodmultiallreduce)(int, N_Vector, sunrealtype*)

      The function implementing :c:func:`N_VDotProdMultiAllReduce`

   .. c:member:: SUNErrCode (*nvbufsize)(N_Vector, sunindextype*)

      The function implementing :c:func:`N_VBufSize`

   .. c:member:: SUNErrCode (*nvbufpack)(N_Vector, void*)

      The function implementing :c:func:`N_VBufPack`

   .. c:member:: SUNErrCode (*nvbufunpack)(N_Vector, void*)

      The function implementing :c:func:`N_VBufUnpack`

   .. c:member:: void (*nvprint)(N_Vector)

      The function implementing :c:func:`N_VPrint`

   .. c:member:: void (*nvprintfile)(N_Vector, FILE*)

      The function implementing :c:func:`N_VPrintFile`

The generic NVECTOR module defines and implements the vector
operations acting on a ``N_Vector``. These routines are nothing but
wrappers for the vector operations defined by a particular NVECTOR
implementation, which are accessed through the *ops* field of the
``N_Vector`` structure. To illustrate this point we show below the
implementation of a typical vector operation from the generic NVECTOR
module, namely ``N_VScale``, which performs the operation :math:`z\gets cx`
for vectors :math:`x` and :math:`z` and a scalar :math:`c`:

.. code-block:: c

   void N_VScale(sunrealtype c, N_Vector x, N_Vector z) {
      z->ops->nvscale(c, x, z);
   }

:numref:`NVectors.Ops` contains a complete list of all standard vector
operations defined by the generic NVECTOR module.  :numref:`NVectors.Ops.Fused`,
:numref:`NVectors.Ops.Array`, :numref:`NVectors.Ops.Local`,
:numref:`NVectors.Ops.SingleBufferReduction`, and
:numref:`NVectors.Ops.Exchange` list *optional* fused, vector array, local
reduction, single buffer reduction, and exchange operations, respectively.

Fused and vector array operations (see :numref:`NVectors.Ops.Fused` and
:numref:`NVectors.Ops.Array`) are intended to increase data reuse, reduce
parallel communication on distributed memory systems, and lower the number of
kernel launches on systems with accelerators. If a particular NVECTOR
implementation defines a fused or vector array operation as ``NULL``, the
generic NVECTOR module will automatically call standard vector operations as
necessary to complete the desired operation. In all SUNDIALS-provided
NVECTOR implementations, all fused and vector array operations are
disabled by default.  However, these implementations provide
additional user-callable functions to enable/disable any or all of the
fused and vector array operations. See the following sections
for the implementation specific functions to enable/disable operations.

Local reduction operations (see :numref:`NVectors.Ops.Local`) are
similarly intended to reduce parallel
communication on distributed memory systems, particularly when
NVECTOR objects are combined together within an NVECTOR_MANYVECTOR
object (see :numref:`NVectors.ManyVector`).  If a
particular NVECTOR implementation defines a local reduction
operation as ``NULL``, the NVECTOR_MANYVECTOR module will
automatically call standard vector reduction operations as necessary
to complete the desired operation. All SUNDIALS-provided NVECTOR
implementations include these local reduction operations, which may be
used as templates for user-defined implementations.

The single buffer reduction operations
(:numref:`NVectors.Ops.SingleBufferReduction`) are used in low-synchronization
methods to combine separate reductions into one ``MPI_Allreduce`` call.

The exchange operations (see :numref:`NVectors.Ops.Exchange`) are intended
only for use with the XBraid library
for parallel-in-time integration (accessible from ARKODE)
and are otherwise unused by SUNDIALS packages.


.. _NVectors.Description.utilities:

NVECTOR Utility Functions
-------------------------

The generic NVECTOR module also defines several utility functions to aid in
creation and management of arrays of ``N_Vector`` objects -- these functions
are particularly useful for Fortran users to utilize the NVECTOR_MANYVECTOR
or SUNDIALS' sensitivity-enabled packages CVODES and IDAS.

The functions :c:func:`N_VCloneVectorArray` and
:c:func:`N_VCloneVectorArrayEmpty` create (by cloning) an array of *count*
variables of type :c:type:`N_Vector`, each of the same type as an existing
``N_Vector`` input:

.. c:function:: N_Vector *N_VCloneVectorArray(int count, N_Vector w)

   Clones an array of ``count``  ``N_Vector`` objects, allocating their data arrays (similar to :c:func:`N_VClone`).

   **Arguments:**
      * ``count`` -- number of ``N_Vector`` objects to create.
      * ``w`` -- template :c:type:`N_Vector` to clone.

   **Return value:**
      * pointer to a new ``N_Vector`` array on success.
      * ``NULL`` pointer on failure.


.. c:function:: N_Vector *N_VCloneVectorArrayEmpty(int count, N_Vector w)

   Clones an array of ``count``  ``N_Vector`` objects, leaving their data arrays unallocated (similar to :c:func:`N_VCloneEmpty`).

   **Arguments:**
      * ``count`` -- number of ``N_Vector`` objects to create.
      * ``w`` -- template :c:type:`N_Vector` to clone.

   **Return value:**
      * pointer to a new ``N_Vector`` array on success.
      * ``NULL`` pointer on failure.


An array of variables of type :c:type:`N_Vector` can be destroyed
by calling :c:func:`N_VDestroyVectorArray`:


.. c:function:: void N_VDestroyVectorArray(N_Vector *vs, int count)

   Destroys an array of ``count``  ``N_Vector`` objects.

   **Arguments:**
      * ``vs`` -- ``N_Vector`` array to destroy.
      * ``count`` -- number of ``N_Vector`` objects in ``vs`` array.

   **Notes:**
      This routine will internally call the ``N_Vector``
      implementation-specific :c:func:`N_VDestroy` operation.

      If ``vs`` was allocated using :c:func:`N_VCloneVectorArray` then
      the data arrays for each ``N_Vector`` object will be freed; if
      ``vs`` was allocated using :c:func:`N_VCloneVectorArrayEmpty` then
      it is the user's responsibility to free the data for each ``N_Vector``
      object.


Finally, we note that users of the Fortran 2003 interface may be interested in
the additional utility functions :c:func:`N_VNewVectorArray`,
:c:func:`N_VGetVecAtIndexVectorArray`, and :c:func:`N_VSetVecAtIndexVectorArray`,
that are wrapped as ``FN_NewVectorArray``, ``FN_VGetVecAtIndexVectorArray``, and
``FN_VSetVecAtIndexVectorArray``, respectively.  These functions allow a Fortran
2003 user to create an empty vector array, access a vector from this array, and
set a vector within this array:


.. c:function:: N_Vector *N_VNewVectorArray(int count, SUNContext sunctx)

   Creates an array of ``count``  ``N_Vector`` objects, the pointers to each
   are initialized as ``NULL``.

   **Arguments:**
      * ``count`` -- length of desired ``N_Vector`` array.
      * ``sunctx`` -- a ``SUNContext`` object

   **Return value:**
      * pointer to a new ``N_Vector`` array on success.
      * ``NULL`` pointer on failure.

   .. versionchanged:: 7.0.0

      The function signature was updated to add the ``SUNContext`` argument.


.. c:function:: N_Vector *N_VGetVecAtIndexVectorArray(N_Vector* vs, int index)

   Accesses the ``N_Vector`` at the location ``index`` within the ``N_Vector`` array ``vs``.

   **Arguments:**
      * ``vs`` -- ``N_Vector`` array.
      * ``index`` -- desired ``N_Vector`` to access from within ``vs``.

   **Return value:**
      * pointer to the indexed ``N_Vector`` on success.
      * ``NULL`` pointer on failure (``index < 0`` or ``vs == NULL``).

   **Notes:**
      This routine does not verify that ``index`` is within the extent of
      ``vs``, since ``vs`` is a simple ``N_Vector`` array that does not
      internally store its allocated length.


.. c:function:: void N_VSetVecAtIndexVectorArray(N_Vector* vs, int index, N_Vector w)

   Sets a pointer to ``w`` at the location ``index`` within the vector array ``vs``.

   **Arguments:**
      * ``vs`` -- ``N_Vector`` array.
      * ``index`` -- desired location to place the pointer to ``w`` within ``vs``.
      * ``w`` -- ``N_Vector`` to set within ``vs``.

   **Notes:**
      This routine does not verify that ``index`` is within the extent of
      ``vs``, since ``vs`` is a simple ``N_Vector`` array that does not
      internally store its allocated length.



.. _NVectors.Description.custom_implementation:

Implementing a custom NVECTOR
-----------------------------

A particular implementation of the NVECTOR module must:

* Specify the *content* field of the ``N_Vector`` structure.

* Define and implement the vector operations.  Note that the names of
  these routines should be unique to that implementation in order to
  permit using more than one NVECTOR module (each with different
  ``N_Vector`` internal data representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free an ``N_Vector`` with
  the new *content* field and with *ops* pointing to the
  new vector operations.

* Optionally, define and implement additional user-callable routines
  acting on the newly-defined ``N_Vector`` (e.g., a routine to print
  the content for debugging purposes).

* Optionally, provide accessor macros as needed for that particular
  implementation to be used to access different parts in the
  *content* field of the newly-defined ``N_Vector``.

To aid in the creation of custom NVECTOR modules, the generic NVECTOR module
provides two utility functions :c:func:`N_VNewEmpty` and
:c:func:`N_VCopyOps()`. When used in custom NVECTOR constructors and clone
routines these functions will ease the introduction of any new optional vector
operations to the NVECTOR API by ensuring that only required operations need
to be set, and that all operations are copied when cloning a vector.

.. c:function:: N_Vector N_VNewEmpty(SUNContext sunctx)

   This allocates a new generic ``N_Vector`` object and initializes its content
   pointer and the function pointers in the operations structure to ``NULL``.

   **Return value:** If successful, this function returns an ``N_Vector``
   object. If an error occurs when allocating the object, then this routine will
   return ``NULL``.

.. c:function:: void N_VFreeEmpty(N_Vector v)

   This routine frees the generic ``N_Vector`` object, under the assumption that any
   implementation-specific data that was allocated within the underlying content structure
   has already been freed. It will additionally test whether the ops pointer is ``NULL``,
   and, if it is not, it will free it as well.

   **Arguments:**
      * *v* -- an N_Vector object

.. c:function:: SUNErrCode N_VCopyOps(N_Vector w, N_Vector v)

   This function copies the function pointers in the ``ops`` structure of ``w``
   into the ``ops`` structure of ``v``.

   **Arguments:**
      * *w* -- the vector to copy operations from
      * *v* -- the vector to copy operations to

   **Return value:**  Returns a :c:type:`SUNErrCode`.


.. c:enum:: N_Vector_ID

   Each :c:type:`N_Vector` implementation included in SUNDIALS has a unique
   identifier specified in enumeration and shown in
   :numref:`NVectors.Description.vectorIDs`. It is recommended that a user
   supplied NVECTOR implementation use the ``SUNDIALS_NVEC_CUSTOM`` identifier.


.. _NVectors.Description.vectorIDs:

.. table:: Vector Identifications associated with vector kernels supplied with SUNDIALS

   ===========================  ====================================  ========
   Vector ID                    Vector type                           ID Value
   ===========================  ====================================  ========
   SUNDIALS_NVEC_SERIAL         Serial                                0
   SUNDIALS_NVEC_PARALLEL       Distributed memory parallel (MPI)     1
   SUNDIALS_NVEC_OPENMP         OpenMP shared memory parallel         2
   SUNDIALS_NVEC_PTHREADS       PThreads shared memory parallel       3
   SUNDIALS_NVEC_PARHYP         *hypre* ParHyp parallel vector        4
   SUNDIALS_NVEC_PETSC          PETSc parallel vector                 5
   SUNDIALS_NVEC_CUDA           CUDA vector                           6
   SUNDIALS_NVEC_HIP            HIP vector                            7
   SUNDIALS_NVEC_SYCL           SYCL vector                           8
   SUNDIALS_NVEC_RAJA           RAJA vector                           9
   SUNDIALS_NVEC_OPENMPDEV      OpenMP vector with device offloading  10
   SUNDIALS_NVEC_TRILINOS       Trilinos Tpetra vector                11
   SUNDIALS_NVEC_MANYVECTOR     "ManyVector" vector                   12
   SUNDIALS_NVEC_MPIMANYVECTOR  MPI-enabled "ManyVector" vector       13
   SUNDIALS_NVEC_MPIPLUSX       MPI+X vector                          14
   SUNDIALS_NVEC_CUSTOM         User-provided custom vector           15
   ===========================  ====================================  ========


.. _NVectors.Description.complex:

Support for complex-valued vectors
----------------------------------

While SUNDIALS itself is written under an assumption of real-valued
data, it does provide limited support for complex-valued problems.
However, since none of the built-in NVECTOR modules supports
complex-valued data, users must provide a custom NVECTOR
implementation for this task.  Many of the NVECTOR routines
described in the subsection :numref:`NVectors.Ops` naturally extend
to complex-valued vectors; however, some do not.  To this end, we
provide the following guidance:

* :c:func:`N_VMin()` and :c:func:`N_VMinLocal()` should return the
  minimum of all *real* components of the vector, i.e.,
  :math:`m = \displaystyle \min_{0\le i< n} \operatorname{real}(x_i)`.

* :c:func:`N_VConst()` (and similarly :c:func:`N_VConstVectorArray()`) should
  set the real components of the vector to the input constant, and set
  all imaginary components to zero, i.e., :math:`z_i = c + 0 j` for :math:`0\le i<n`.

* :c:func:`N_VAddConst()` should only update the real components of the
  vector with the input constant, leaving all imaginary components
  unchanged.

* :c:func:`N_VWrmsNorm()`, :c:func:`N_VWrmsNormMask()`,
  :c:func:`N_VWSqrSumLocal()` and :c:func:`N_VWSqrSumMaskLocal()`
  should assume that all entries of the weight vector ``w`` and the
  mask vector ``id`` are real-valued.

* :c:func:`N_VDotProd()` should mathematically return a complex number
  for complex-valued vectors; as this is not possible with
  SUNDIALS' current ``sunrealtype``, this routine should
  be set to ``NULL`` in the custom NVECTOR implementation.

* :c:func:`N_VCompare()`, :c:func:`N_VConstrMask()`, :c:func:`N_VMinQuotient()`,
  :c:func:`N_VConstrMaskLocal()` and :c:func:`N_VMinQuotientLocal()`
  are ill-defined due to the lack of a clear ordering in the
  complex plane.  These routines should be set to ``NULL``
  in the custom NVECTOR implementation.


While many SUNDIALS solver modules may be utilized on complex-valued data,
others cannot.  Specifically, although each package's linear solver
interface (e.g., ARKLS or CVLS) may be used on complex-valued problems,
none of the built-in SUNMatrix or SUNLinearSolver modules will work (all
of the direct linear solvers must store complex-valued data, and all of
the iterative linear solvers require :c:func:`N_VDotProd`).  Hence a
complex-valued user must provide custom linear solver modules for their
problem.  At a minimum this will consist of a custom SUNLinearSolver
implementation (see :numref:`SUNLinSol.API.Custom`), and optionally a
custom SUNMatrix as well.  The user should then attach these modules as
normal to the package's linear solver interface.

.. ifconfig:: package_name != 'kinsol'

   Similarly, although both the
   :ref:`SUNNonlinearSolver_Newton <SUNNonlinSol.Newton>` and
   :ref:`SUNNonlinearSolver_FixedPoint <SUNNonlinSol.FixedPoint>` modules
   may be used with any of the IVP solvers (CVODE(S), IDA(S) and ARKODE) for
   complex-valued problems, the Anderson-acceleration option with
   SUNNonlinearSolver_FixedPoint cannot be used due to its reliance on
   :c:func:`N_VDotProd()`.  By this same logic, the Anderson acceleration
   feature within KINSOL will also not work with complex-valued vectors.

Finally, constraint-handling features of each package cannot be used
for complex-valued data, due to the issue of
ordering in the complex plane discussed above with
:c:func:`N_VCompare()`, :c:func:`N_VConstrMask()`,
:c:func:`N_VMinQuotient()`, :c:func:`N_VConstrMaskLocal()` and
:c:func:`N_VMinQuotientLocal()`.

We provide a simple example of a complex-valued example problem,
including a custom complex-valued Fortran 2003 NVECTOR module, in the
files ``examples/arkode/F2003_custom/ark_analytic_complex_f2003.f90``,
``examples/arkode/F2003_custom/fnvector_complex_mod.f90``, and
``examples/arkode/F2003_custom/test_fnvector_complex_mod.f90``.
