! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2025-2026, Lawrence Livermore National Security,
! University of Maryland Baltimore County, and the SUNDIALS contributors.
! Copyright (c) 2013-2025, Lawrence Livermore National Security
! and Southern Methodist University.
! Copyright (c) 2002-2013, Lawrence Livermore National Security.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------
! This file tests the Fortran 2003 interface to the SUNDIALS
! manyvector N_Vector implementation.
! -----------------------------------------------------------------

module test_nvector_manyvector
  use, intrinsic :: iso_c_binding

  use fnvector_manyvector_mod
  use fnvector_serial_mod
  use test_utilities
  implicit none

  integer(c_int), parameter  :: nsubvecs = 2
  integer(c_int), parameter  :: nv = 3       ! length of vector arrays
  integer(kind=myindextype), parameter  :: N1 = 100     ! individual vector length
  integer(kind=myindextype), parameter  :: N2 = 200     ! individual vector length
  integer(kind=myindextype), parameter  :: N = N1 + N2 ! overall manyvector length

contains

  integer function smoke_tests() result(ret)
    implicit none

    integer(kind=myindextype) :: ival                ! integer work value
    integer(kind=myindextype) :: lenrw(1), leniw(1)  ! real and int work space size
    real(c_double)          :: rval                   ! real work value
    real(c_double)          :: x1data(N1), x2data(N2) ! vector data array
    real(c_double), pointer :: xptr(:)                ! pointer to vector data array
    real(c_double)          :: nvarr(nv)              ! array of nv constants to go with vector array
    type(N_Vector), pointer :: x, y, z, tmp           ! N_Vectors
    type(c_ptr)             :: subvecs                ! ManyVector subvectors
    type(c_ptr)             :: xvecs, zvecs           ! C pointer to array of ManyVectors

    !===== Setup ====
    subvecs = FN_VNewVectorArray(nsubvecs, sunctx)
    tmp => FN_VMake_Serial(N1, x1data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 0, tmp)
    tmp => FN_VMake_Serial(N2, x2data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 1, tmp)

    x => FN_VNew_ManyVector(int(nsubvecs, myindextype), subvecs, sunctx)
    call FN_VConst(ONE, x)
    y => FN_VClone_ManyVector(x)
    call FN_VConst(ONE, y)
    z => FN_VClone_ManyVector(x)
    call FN_VConst(ONE, z)

    xvecs = FN_VCloneVectorArray(nv, x)
    zvecs = FN_VCloneVectorArray(nv, z)
    nvarr = (/ONE, ONE, ONE/)

    !===== Test =====

    ! test generic vector functions
    ival = FN_VGetVectorID_ManyVector(x)
    call FN_VSpace_ManyVector(x, lenrw, leniw)
    ival = FN_VGetCommunicator(x)
    ival = FN_VGetLength_ManyVector(x)

    ! test standard vector operations
    call FN_VLinearSum_ManyVector(ONE, x, ONE, y, z)
    call FN_VConst_ManyVector(ONE, z)
    call FN_VProd_ManyVector(x, y, z)
    call FN_VDiv_ManyVector(x, y, z)
    call FN_VScale_ManyVector(ONE, x, y)
    call FN_VAbs_ManyVector(x, y)
    call FN_VInv_ManyVector(x, z)
    call FN_VAddConst_ManyVector(x, ONE, z)
    rval = FN_VDotProdLocal_ManyVector(x, y)
    rval = FN_VMaxNormLocal_ManyVector(x)
    rval = FN_VWrmsNorm_ManyVector(x, y)
    rval = FN_VWrmsNormMask_ManyVector(x, y, z)
    rval = FN_VMinLocal_ManyVector(x)
    rval = FN_VWL2Norm_ManyVector(x, y)
    rval = FN_VL1NormLocal_ManyVector(x)
    call FN_VCompare_ManyVector(ONE, x, y)
    ival = FN_VInvTestLocal_ManyVector(x, y)
    ival = FN_VConstrMaskLocal_ManyVector(z, x, y)
    rval = FN_VMinQuotientLocal_ManyVector(x, y)

    ! test fused vector operations
    ival = FN_VLinearCombination_ManyVector(nv, nvarr, xvecs, x)
    ival = FN_VScaleAddMulti_ManyVector(nv, nvarr, x, xvecs, zvecs)
    ival = FN_VDotProdMulti_ManyVector(nv, x, xvecs, nvarr)

    ! test vector array operations
    ival = FN_VLinearSumVectorArray_ManyVector(nv, ONE, xvecs, ONE, xvecs, zvecs)
    ival = FN_VScaleVectorArray_ManyVector(nv, nvarr, xvecs, zvecs)
    ival = FN_VConstVectorArray_ManyVector(nv, ONE, xvecs)
    ival = FN_VWrmsNormVectorArray_ManyVector(nv, xvecs, xvecs, nvarr)
    ival = FN_VWrmsNormMaskVectorArray_ManyVector(nv, xvecs, xvecs, x, nvarr)

    ! test the ManyVector specific operations
    ival = FN_VGetNumSubvectors_ManyVector(x)
    xptr => FN_VGetSubvectorArrayPointer_ManyVector(x, ival - 1)
    ival = FN_VSetSubvectorArrayPointer_ManyVector(xptr, x, ival - 1)
    ival = FN_VGetNumSubvectors_ManyVector(x)
    tmp => FN_VGetSubvector_ManyVector(x, ival - 1)

    !==== Cleanup =====
    call FN_VDestroyVectorArray(subvecs, nsubvecs)
    call FN_VDestroy_ManyVector(x)
    call FN_VDestroy_ManyVector(y)
    call FN_VDestroy_ManyVector(z)
    call FN_VDestroyVectorArray(xvecs, nv)
    call FN_VDestroyVectorArray(zvecs, nv)

    ret = 0

  end function smoke_tests

  integer function unit_tests() result(fails)
    use test_fnvector
    implicit none

    real(c_double)           :: x1data(N1)   ! vector data array
    real(c_double)           :: x2data(N2)   ! vector data array
    type(N_Vector), pointer  :: x, tmp       ! N_Vectors
    type(c_ptr)              :: subvecs      ! subvectors of the ManyVector

    !===== Setup ====
    fails = 0

    subvecs = FN_VNewVectorArray(nsubvecs, sunctx)
    tmp => FN_VMake_Serial(N1, x1data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 0, tmp)
    tmp => FN_VMake_Serial(N2, x2data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 1, tmp)

    x => FN_VNew_ManyVector(int(nsubvecs, myindextype), subvecs, sunctx)
    call FN_VConst(ONE, x)

    !==== tests ====
    fails = Test_FN_VMake(x, N, 0)
    fails = Test_FN_VGetArrayPointer(x, N, 0)
    fails = Test_FN_VLinearCombination(x, N, 0)

    !=== cleanup ====
    call FN_VDestroyVectorArray(subvecs, nsubvecs)
    call FN_VDestroy_ManyVector(x)

  end function unit_tests

end module

function check_ans(ans, X, local_length) result(failure)
  use, intrinsic :: iso_c_binding
  use fnvector_manyvector_mod
  use test_utilities
  implicit none

  real(C_DOUBLE)             :: ans
  type(N_Vector)             :: X
  type(N_Vector), pointer    :: X0, X1
  integer(kind=myindextype) :: failure, local_length, i, x0len, x1len
  real(C_DOUBLE), pointer    :: x0data(:), x1data(:)

  failure = 0

  X0 => FN_VGetSubvector_ManyVector(X, 0_myindextype)
  X1 => FN_VGetSubvector_ManyVector(X, 1_myindextype)
  x0len = FN_VGetLength(X0)
  x1len = FN_VGetLength(X1)
  x0data => FN_VGetArrayPointer(X0)
  x1data => FN_VGetArrayPointer(X1)

  if (local_length /= (x0len + x1len)) then
    failure = 1
    return
  end if

  do i = 1, x0len
    if (FNEQ(x0data(i), ans) > 0) then
      failure = failure + 1
    end if
  end do

  do i = 1, x1len
    if (FNEQ(x1data(i), ans) > 0) then
      failure = failure + 1
    end if
  end do
end function check_ans

logical function has_data(X) result(failure)
  use, intrinsic :: iso_c_binding

  use test_utilities
  implicit none

  type(N_Vector)          :: X

  failure = .true.
end function has_data

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_nvector_manyvector

  !======== Declarations ========
  implicit none
  integer :: fails = 0

  !============== Introduction =============
  print *, 'ManyVector N_Vector Fortran 2003 interface test'

  call Test_Init(SUN_COMM_NULL)

  fails = smoke_tests()
  if (fails /= 0) then
    print *, 'FAILURE: smoke tests failed'
    stop 1
  else
    print *, 'SUCCESS: smoke tests passed'
  end if

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed'
    stop 1
  else
    print *, 'SUCCESS: all unit tests passed'
  end if

  call Test_Finalize()
end program main
