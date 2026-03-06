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

%module fnvector_manyvector_mod

// include code common to all nvector implementations
%include "fnvector.i"

// include the header file in the swig wrapper
%{
#include "nvector/nvector_manyvector.h"
%}

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(ManyVector)

// Process and wrap functions in the following files
%include "nvector/nvector_manyvector.h"

%insert("wrapper") %{
#ifdef SUNDIALS_INT32_T
SWIGEXPORT double * _wrap_FN_VGetSubvectorArrayPointer_ManyVector(N_Vector farg1, int32_t const *farg2) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunindextype arg2 ;
  sunrealtype *result = 0 ;

  arg1 = (N_Vector)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (sunrealtype *)N_VGetSubvectorArrayPointer_ManyVector(arg1,arg2);
  fresult = result;
  return fresult;
}
#else
SWIGEXPORT double * _wrap_FN_VGetSubvectorArrayPointer_ManyVector(N_Vector farg1, int64_t const *farg2) {
  double * fresult ;
  N_Vector arg1 = (N_Vector) 0 ;
  sunindextype arg2 ;
  sunrealtype *result = 0 ;

  arg1 = (N_Vector)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (sunrealtype *)N_VGetSubvectorArrayPointer_ManyVector(arg1,arg2);
  fresult = result;
  return fresult;
}
#endif
%}

%insert("fdecl") %{
 public :: FN_VGetSubvectorArrayPointer_ManyVector
%}

%insert("finterfaces") %{
function swigc_FN_VGetSubvectorArrayPointer_ManyVector(farg1, farg2) &
bind(C, name="_wrap_FN_VGetSubvectorArrayPointer_ManyVector") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T), intent(in) :: farg2
#else
integer(C_INT64_T), intent(in) :: farg2
#endif
type(C_PTR) :: fresult
end function
%}

%insert("fsubprograms") %{
function FN_VGetSubvectorArrayPointer_ManyVector(v, vec_num) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(N_Vector), target, intent(inout) :: v
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T), intent(in) :: vec_num
#else
integer(C_INT64_T), intent(in) :: vec_num
#endif
type(C_PTR) :: fresult
type(C_PTR) :: farg1
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T) :: farg2
#else
integer(C_INT64_T) :: farg2
#endif

farg1 = c_loc(v)
farg2 = vec_num
fresult = swigc_FN_VGetSubvectorArrayPointer_ManyVector(farg1, farg2)
call c_f_pointer(fresult, swig_result, [FN_VGetSubvectorLocalLength_ManyVector(v, vec_num)])
end function
%}
