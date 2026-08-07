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

module test_nvector_mpiplusx
  use, intrinsic :: iso_c_binding

  use fnvector_mpiplusx_mod
  use fnvector_serial_mod
  use test_utilities
  implicit none
  include "mpif.h"

  integer(kind=myindextype), parameter :: N = 100     ! overall manyvector length
  integer(c_int), target     :: comm = MPI_COMM_WORLD ! default MPI communicator
  integer(c_int)             :: nprocs                ! number of MPI processes

contains

  integer function smoke_tests() result(ret)
    implicit none

    integer(c_long)         :: ival           ! integer work value
    real(c_double)          :: x1data(N)     ! vector data array
    real(c_double), pointer :: xptr(:)        ! pointer to vector data array
    type(N_Vector), pointer :: x, local       ! N_Vectors

    !===== Setup ====
    local => FN_VMake_Serial(N, x1data, sunctx)

    x => FN_VMake_MPIPlusX(comm, local, sunctx)
    call FN_VConst(ONE, x)

    !===== Test =====

    ! test the MPIPlusX specific operations
    xptr => FN_VGetArrayPointer_MPIPlusX(x)
    local => FN_VGetLocalVector_MPIPlusX(x)
    ival = FN_VGetLocalLength_MPIPlusX(x)
    ival = FN_VGetVectorID_MPIPlusX(x)

    !==== Cleanup =====
    call FN_VDestroy(local)
    call FN_VDestroy(x)

    ret = 0

  end function smoke_tests

  integer function unit_tests() result(fails)
    use test_fnvector
    implicit none

    integer(c_int)           :: myid         ! the MPI rank
    real(c_double)           :: x1data(N)   ! vector data array
    type(N_Vector), pointer  :: x, local       ! N_Vectors

    !===== Setup ====
    fails = 0

    call MPI_Comm_rank(comm, myid, fails)
    if (fails /= 0) then
      print *, '   FAILURE - MPI_COMM_RANK returned nonzero'
      stop 1
    end if

    local => FN_VMake_Serial(N, x1data, sunctx)
    x => FN_VMake_MPIPlusX(comm, local, sunctx)
    call FN_VConst(ONE, x)

    !==== tests ====
    fails = Test_FN_VMake(x, N, myid)
    fails = Test_FN_VGetArrayPointer(x, N, myid)
    fails = Test_FN_VLinearCombination(x, N, myid)

    !=== cleanup ====
    local => FN_VGetLocalVector_MPIPlusX(x)
    call FN_VDestroy(local)
    call FN_VDestroy(x)

  end function unit_tests

end module

integer(C_INT) function check_ans(ans, X, local_length) result(failure)
  use, intrinsic :: iso_c_binding
  use fnvector_mpiplusx_mod

  use test_utilities
  implicit none

  real(C_DOUBLE)            :: ans
  type(N_Vector)            :: X
  type(N_Vector), pointer   :: X0
  integer(kind=myindextype) :: local_length, i, x0len
  real(C_DOUBLE), pointer   :: x0data(:)

  failure = 0

  X0 => FN_VGetLocalVector_MPIPlusX(X)
  x0len = FN_VGetLength(X0)
  x0data => FN_VGetArrayPointer(X0)

  if (local_length /= x0len) then
    failure = 1
    return
  end if

  do i = 1, x0len
    if (FNEQ(x0data(i), ans) > 0) then
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
  use test_nvector_mpiplusx

  !======== Declarations ========
  implicit none
  integer(c_int) :: fails = 0
  integer(c_int) :: myid

  !===== MPI Setup ====
  call MPI_Init(fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_INIT returned nonzero'
    stop 1
  end if

  call MPI_Comm_rank(comm, myid, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_COMM_RANK returned nonzero, proc', myid
    stop 1
  end if

  !============== Introduction =============
  if (myid == 0) print *, 'MPIPlusX N_Vector Fortran 2003 interface test'

  call Test_Init(comm)

  call MPI_Comm_size(comm, nprocs, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_COMM_SIZE returned nonzero, proc', myid
    stop 1
  end if

  fails = smoke_tests()
  if (fails /= 0) then
    print *, 'FAILURE: smoke tests failed, proc', myid
  else
    if (myid == 0) &
      print *, 'SUCCESS: all smoke tests passed'
  end if

  call MPI_BARRIER(comm, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_BARRIER returned nonzero, proc', myid
    stop 1
  end if

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed, proc', myid
    stop 1
  else
    if (myid == 0) print *, '    SUCCESS - all unit tests passed'
  end if

  call Test_Finalize()

  call MPI_Finalize(fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_FINALIZE returned nonzero, proc ', myid
    stop 1
  end if
end program main
