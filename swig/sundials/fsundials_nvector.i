// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2025-2026, Lawrence Livermore National Security,
// University of Maryland Baltimore County, and the SUNDIALS contributors.
// Copyright (c) 2013-2025, Lawrence Livermore National Security
// and Southern Methodist University.
// Copyright (c) 2002-2013, Lawrence Livermore National Security.
// All rights reserved.
//
// See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-3-Clause
// SUNDIALS Copyright End
// ---------------------------------------------------------------
// Swig interface file
// ---------------------------------------------------------------


// insert the include into the swig wrapper
%{
#include "sundials/sundials_nvector.h"
%}

// Ignore functions with arrays of vector arrays since they are not supported
%ignore N_VScaleAddMultiVectorArray;
%ignore N_VLinearCombinationVectorArray;

// Ignore GetArrayPointer functions because we manually insert them
%ignore N_VGetArrayPointer;
%ignore N_VGetDeviceArrayPointer;

// Process and wrap functions in the following files
%include "sundials/sundials_nvector.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FN_VGetArrayPointer(N_Vector farg1) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunrealtype *result = 0 ;
  arg1 = (N_Vector)(farg1);
  result = (sunrealtype *)N_VGetArrayPointer(arg1);
  fresult = result;
  return fresult;
}


SWIGEXPORT double * _wrap_FN_VGetDeviceArrayPointer(N_Vector farg1) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunrealtype *result = 0 ;
  arg1 = (N_Vector)(farg1);
  result = (sunrealtype *)N_VGetDeviceArrayPointer(arg1);
  fresult = result;
  return fresult;
}
%}

%insert("fdecl") %{
 public :: FN_VGetArrayPointer
 public :: FN_VGetDeviceArrayPointer
%}

%insert("finterfaces") %{
function swigc_FN_VGetArrayPointer(farg1) &
bind(C, name="_wrap_FN_VGetArrayPointer") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FN_VGetDeviceArrayPointer(farg1) &
bind(C, name="_wrap_FN_VGetDeviceArrayPointer") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function
%}

%insert("fsubprograms") %{
function FN_VGetArrayPointer(v) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(N_Vector), target, intent(inout) :: v
type(C_PTR) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(v)
fresult = swigc_FN_VGetArrayPointer(farg1)
call c_f_pointer(fresult, swig_result, [FN_VGetLocalLength(v)])
end function

function FN_VGetDeviceArrayPointer(v) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(N_Vector), target, intent(inout) :: v
type(C_PTR) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(v)
fresult = swigc_FN_VGetDeviceArrayPointer(farg1)
call c_f_pointer(fresult, swig_result, [FN_VGetLocalLength(v)])
end function
%}
