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

%module fsunmatrix_sparse_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

%{
#include "sunmatrix/sunmatrix_sparse.h"
%}

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Sparse)

// we ignore the following functions because we manually insert them to ensure the
// the returned arrays have the correct shape
%ignore SUNSparseMatrix_Data;
%ignore SUNSparseMatrix_IndexValues;
%ignore SUNSparseMatrix_IndexPointers;

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_sparse.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FSUNSparseMatrix_Data(SUNMatrix farg1) {
  double * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunrealtype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  result = (sunrealtype *)SUNSparseMatrix_Data(arg1);
  fresult = result;
  return fresult;
}

#ifdef SUNDIALS_INT32_T
SWIGEXPORT int32_t * _wrap_FSUNSparseMatrix_IndexValues(SUNMatrix farg1) {
  int32_t * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  result = (sunindextype *)SUNSparseMatrix_IndexValues(arg1);
  fresult = result;
  return fresult;
}

SWIGEXPORT int32_t * _wrap_FSUNSparseMatrix_IndexPointers(SUNMatrix farg1) {
  int32_t * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  result = (sunindextype *)SUNSparseMatrix_IndexPointers(arg1);
  fresult = result;
  return fresult;
}
#else
SWIGEXPORT int64_t * _wrap_FSUNSparseMatrix_IndexValues(SUNMatrix farg1) {
  int64_t * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  result = (sunindextype *)SUNSparseMatrix_IndexValues(arg1);
  fresult = result;
  return fresult;
}

SWIGEXPORT int64_t * _wrap_FSUNSparseMatrix_IndexPointers(SUNMatrix farg1) {
  int64_t * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  result = (sunindextype *)SUNSparseMatrix_IndexPointers(arg1);
  fresult = result;
  return fresult;
}
#endif
%}

%insert("fdecl") %{
 public :: FSUNSparseMatrix_Data
 public :: FSUNSparseMatrix_IndexValues
 public :: FSUNSparseMatrix_IndexPointers
%}

%insert("finterfaces") %{
function swigc_FSUNSparseMatrix_Data(farg1) &
bind(C, name="_wrap_FSUNSparseMatrix_Data") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FSUNSparseMatrix_IndexValues(farg1) &
bind(C, name="_wrap_FSUNSparseMatrix_IndexValues") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FSUNSparseMatrix_IndexPointers(farg1) &
bind(C, name="_wrap_FSUNSparseMatrix_IndexPointers") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function
%}

%insert("fsubprograms") %{

function FSUNSparseMatrix_Data(a) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(SUNMatrix), target, intent(inout) :: a
type(C_PTR) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(a)
fresult = swigc_FSUNSparseMatrix_Data(farg1)
call c_f_pointer(fresult, swig_result, [FSUNSparseMatrix_NNZ(a)])
end function

function FSUNSparseMatrix_IndexValues(a) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T), dimension(:), pointer :: swig_result
#else
integer(C_INT64_T), dimension(:), pointer :: swig_result
#endif
type(SUNMatrix), target, intent(inout) :: a
type(C_PTR) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(a)
fresult = swigc_FSUNSparseMatrix_IndexValues(farg1)
call c_f_pointer(fresult, swig_result, [FSUNSparseMatrix_NNZ(a)])
end function

function FSUNSparseMatrix_IndexPointers(a) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T), dimension(:), pointer :: swig_result
#else
integer(C_INT64_T), dimension(:), pointer :: swig_result
#endif
type(SUNMatrix), target, intent(inout) :: a
type(C_PTR) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(a)
fresult = swigc_FSUNSparseMatrix_IndexPointers(farg1)
call c_f_pointer(fresult, swig_result, [FSUNSparseMatrix_NP(a)+1])
end function
%}
