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

%module fsunmatrix_band_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

%{
#include "sunmatrix/sunmatrix_band.h"
%}

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Band)

// we manually insert these so that the correct shape array is returned
%ignore SUNBandMatrix_Data;
%ignore SUNBandMatrix_Column;

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_band.h"

%insert("wrapper") %{
SWIGEXPORT double * _wrap_FSUNBandMatrix_Data(SUNMatrix farg1) {
  double * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunrealtype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  result = (sunrealtype *)SUNBandMatrix_Data(arg1);
  fresult = result;
  return fresult;
}

#ifdef SUNDIALS_INT32_T
SWIGEXPORT double * _wrap_FSUNBandMatrix_Column(SUNMatrix farg1, int32_t const *farg2) {
  double * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype arg2 ;
  sunrealtype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (sunrealtype *)SUNBandMatrix_Column(arg1,arg2);
  fresult = result;
  return fresult;
}
#else
SWIGEXPORT double * _wrap_FSUNBandMatrix_Column(SUNMatrix farg1, int64_t const *farg2) {
  double * fresult ;
  SUNMatrix arg1 = (SUNMatrix) 0 ;
  sunindextype arg2 ;
  sunrealtype *result = 0 ;

  arg1 = (SUNMatrix)(farg1);
  arg2 = (sunindextype)(*farg2);
  result = (sunrealtype *)SUNBandMatrix_Column(arg1,arg2);
  fresult = result;
  return fresult;
}
#endif
%}

%insert("fdecl") %{
 public :: FSUNBandMatrix_Data
 public :: FSUNBandMatrix_Column
%}

%insert("finterfaces") %{
function swigc_FSUNBandMatrix_Data(farg1) &
bind(C, name="_wrap_FSUNBandMatrix_Data") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FSUNBandMatrix_Column(farg1, farg2) &
bind(C, name="_wrap_FSUNBandMatrix_Column") &
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
function FSUNBandMatrix_Data(a) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(SUNMatrix), target, intent(inout) :: a
type(C_PTR) :: fresult
type(C_PTR) :: farg1

farg1 = c_loc(a)
fresult = swigc_FSUNBandMatrix_Data(farg1)
call c_f_pointer(fresult, swig_result, [FSUNBandMatrix_LData(a)])
end function

function FSUNBandMatrix_Column(a, j) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(:), pointer :: swig_result
type(SUNMatrix), target, intent(inout) :: a
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T), intent(in) :: j
#else
integer(C_INT64_T), intent(in) :: j
#endif
type(C_PTR) :: fresult
type(C_PTR) :: farg1
#ifdef SUNDIALS_INT32_T
integer(C_INT32_T) :: farg2
#else
integer(C_INT64_T) :: farg2
#endif

farg1 = c_loc(a)
farg2 = j
fresult = swigc_FSUNBandMatrix_Column(farg1, farg2)
! We set the array shape to [1] because only the diagonal element
! can be accessed through this function from Fortran.
call c_f_pointer(fresult, swig_result, [1])
end function
%}

