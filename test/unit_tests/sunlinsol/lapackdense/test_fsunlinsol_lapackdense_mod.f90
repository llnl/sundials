! -----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ UMBC
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
! LAPACK-Dense SUNLinearSolver implementation.
! -----------------------------------------------------------------

module test_fsunlinsol_lapackdense
  use, intrinsic :: iso_c_binding
  use test_utilities
  implicit none

  integer(kind=myindextype), private, parameter :: N = 100

contains

  integer(c_int) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding

    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsunlinsol_lapackdense_mod
    use test_sunlinsol

    implicit none

    type(SUNLinearSolver), pointer :: LS                ! test linear solver
    type(SUNMatrix), pointer :: A, I              ! test matrices
    type(N_Vector), pointer :: x, b              ! test vectors
    real(c_double), pointer :: colj(:), colIj(:) ! matrix column data
    real(c_double), pointer :: xdata(:)          ! x vector data
    real(c_double)                 :: tmpr              ! temporary real value
    integer(kind=myindextype)     :: j, k
    integer(c_int)                 :: tmp

    fails = 0

    A => FSUNDenseMatrix(N, N, sunctx)
    I => FSUNDenseMatrix(N, N, sunctx)
    x => FN_VNew_Serial(N, sunctx)
    b => FN_VNew_Serial(N, sunctx)

    ! fill A matrix with uniform random data in [0, 1/N)
    do j = 1, N
      colj => FSUNDenseMatrix_Column(A, j - 1)
      do k = 1, N
        call random_number(tmpr)
        colj(k) = tmpr/N
      end do
    end do

    ! create anti-identity matrix
    j = N
    do k = 1, N
      colj => FSUNDenseMatrix_Column(I, j - 1)
      colj(k) = ONE
      j = j - 1
    end do

    ! add anti-identity to ensure the solver needs to do row-swapping
    do k = 1, N
      do j = 1, N
        colj => FSUNDenseMatrix_Column(A, j - 1)
        colIj => FSUNDenseMatrix_Column(I, j - 1)
        colj(k) = colj(k) + colIj(k)
      end do
    end do

    ! fill x vector with uniform random data in [0, 1)
    xdata => FN_VGetArrayPointer(x)
    do j = 1, N
      call random_number(tmpr)
      xdata(j) = tmpr
    end do

    ! create RHS vector for linear solve
    fails = FSUNMatMatvec(A, x, b)
    if (fails /= 0) then
      call FSUNMatDestroy(A)
      call FSUNMatDestroy(I)
      call FN_VDestroy(x)
      call FN_VDestroy(b)
      return
    end if

    ! create dense linear solver
    LS => FSUNLinSol_LapackDense(x, A, sunctx)

    ! run tests
    fails = fails + Test_FSUNLinSolInitialize(LS, 0)
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0)
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, 100*SUN_UNIT_ROUNDOFF, 0)

    fails = fails + Test_FSUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0)
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0)
    fails = fails + Test_FSUNLinSolSpace(LS, 0)

    ! cleanup
    tmp = FSUNLinSolFree(LS)
    call FSUNMatDestroy(A)
    call FSUNMatDestroy(I)
    call FN_VDestroy(x)
    call FN_VDestroy(b)

  end function unit_tests

end module

integer(c_int) function check_vector(X, Y, tol) result(failure)
  use, intrinsic :: iso_c_binding
  use test_utilities

  implicit none
  type(N_Vector)  :: x, y
  real(c_double)  :: tol, maxerr
  integer(c_long) :: i, xlen, ylen
  real(c_double), pointer :: xdata(:), ydata(:)

  failure = 0

  xdata => FN_VGetArrayPointer(x)
  ydata => FN_VGetArrayPointer(y)

  xlen = FN_VGetLength(x)
  ylen = FN_VGetLength(y)

  if (xlen /= ylen) then
    print *, 'FAIL: check_vector: different data array lengths'
    failure = 1
    return
  end if

  do i = 1, xlen
    failure = failure + FNEQTOL(xdata(i), ydata(i), tol)
  end do

  if (failure > 0) then
    maxerr = ZERO
    do i = 1, xlen
      maxerr = max(abs(xdata(i) - ydata(i))/abs(ydata(i)), maxerr)
    end do
    write (*, '(A,E14.7,A,E14.7,A)') &
      "FAIL: check_vector failure: maxerr = ", maxerr, "  (tol = ", tol, ")"
  end if

end function check_vector

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunlinsol_lapackdense

  !======== Declarations ========
  implicit none
  integer(c_int) :: fails = 0

  !============== Introduction =============
  print *, 'LAPACK-Dense SUNLinearSolver Fortran 2003 interface test'

  call Test_Init(SUN_COMM_NULL)

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed'
    stop 1
  else
    print *, 'SUCCESS: all unit tests passed'
  end if

  call Test_Finalize()

end program main
