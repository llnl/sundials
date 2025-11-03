.. ----------------------------------------------------------------
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

.. _Fortran.DataTypes:

Data Types
==========

Generally, the Fortran 2003 type that is equivalent to the C type is what one
would expect. Primitive types map to the ``iso_c_binding`` type equivalent.
SUNDIALS classes map to a Fortran derived type. However, the handling of pointer
types is not always clear as they can depend on the parameter direction.
:numref:`Fortran.DataTypes.Table` presents a summary of the type
equivalencies with the parameter direction in mind.

.. warning::

   Currently, the Fortran 2003 interfaces are only compatible with SUNDIALS
   builds where the ``sunrealtype`` is double-precision.

.. versionchanged:: 7.1.0

   The Fortran interfaces can now be built with 32-bit ``sunindextype`` in
   addition to 64-bit ``sunindextype``.


.. _Fortran.DataTypes.Table:
.. list-table:: C/Fortran-2003 Equivalent Types. ``T`` represents any type.
   :align: center
   :header-rows: 1

   * - **C Type**
     - **Parameter Direction**
     - **Fortran 2003 type**
   * - ``SUNComm``
     - in, inout, out, return
     - ``integer(c_int)``
   * - ``SUNErrCode``
     - in, inout, out, return
     - ``integer(c_int)``
   * - ``double``
     - in, inout, out, return
     - ``real(c_double)``
   * - ``int``
     - in, inout, out, return
     - ``integer(c_int)``
   * - ``long``
     - in, inout, out, return
     - ``integer(c_long)``
   * - ``sunbooleantype``
     - in, inout, out, return
     - ``integer(c_int)``
   * - ``sunrealtype``
     - in, inout, out, return
     - ``real(c_double)``
   * - ``sunindextype``
     - in, inout, out, return
     - ``integer(c_long)``
   * - ``double*``
     - in, inout, out
     - ``real(c_double), dimension(*)``
   * - ``double*``
     - return
     - ``real(c_double), pointer, dimension(:)``
   * - ``int*``
     - in, inout, out
     - ``real(c_int), dimension(*)``
   * - ``int*``
     - return
     - ``real(c_int), pointer, dimension(:)``
   * - ``long*``
     - in, inout, out
     - ``real(c_long), dimension(*)``
   * - ``long*``
     - return
     - ``real(c_long), pointer, dimension(:)``
   * - ``sunrealtype*``
     - in, inout, out
     - ``real(c_double), dimension(*)``
   * - ``sunrealtype*``
     - return
     - ``real(c_double), pointer, dimension(:)``
   * - ``sunindextype*``
     - in, inout, out
     - ``real(c_long), dimension(*)``
   * - ``sunindextype*``
     - return
     - ``real(c_long), pointer, dimension(:)``
   * - ``sunrealtype[]``
     - in, inout, out
     - ``real(c_double), dimension(*)``
   * - ``sunindextype[]``
     - in, inout, out
     - ``integer(c_long), dimension(*)``
   * - ``SUNAdaptController``
     - in, inout, out
     - ``type(SUNAdaptController)``
   * - ``SUNAdaptController``
     - return
     - ``type(SUNAdaptController), pointer``
   * - ``SUNAdjointCheckpointScheme``
     - in, inout, out, return
     - ``type(c_ptr)``
   * - ``SUNAdjointStepper``
     - in, inout, out, return
     - ``type(c_ptr)``
   * - ``SUNDomEigEstimator``
     - in, inout, out
     - ``type(SUNDomEigEstimator)``
   * - ``SUNDomEigEstimator``
     - return
     - ``type(SUNDomEigEstimator), pointer``
   * - ``SUNMatrix``
     - in, inout, out
     - ``type(SUNMatrix)``
   * - ``SUNMatrix``
     - return
     - ``type(SUNMatrix), pointer``
   * - ``SUNLinearSolver``
     - in, inout, out
     - ``type(SUNLinearSolver)``
   * - ``SUNLinearSolver``
     - return
     - ``type(SUNLinearSolver), pointer``
   * - ``SUNNonlinearSolver``
     - in, inout, out
     - ``type(SUNNonlinearSolver)``
   * - ``SUNNonlinearSolver``
     - return
     - ``type(SUNNonlinearSolver), pointer``
   * - ``N_Vector``
     - in, inout, out
     - ``type(N_Vector)``
   * - ``N_Vector``
     - return
     - ``type(N_Vector), pointer``
   * - ``FILE*``
     - in, inout, out, return
     - ``type(c_ptr)``
   * - ``void*``
     - in, inout, out, return
     - ``type(c_ptr)``
   * - ``T**``
     - in, inout, out, return
     - ``type(c_ptr)``
   * - ``T***``
     - in, inout, out, return
     - ``type(c_ptr)``
   * - ``T****``
     - in, inout, out, return
     - ``type(c_ptr)``
