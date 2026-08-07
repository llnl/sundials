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
! SPGMR SUNLinearSolver implementation. We run three tests to
! exercise the solver:
!   1. simple tridiagonal system (no preconditioning)
!   2. simple tridiagonal system (jacobi preconditioning)
!   3. tridiagonal system w/ scale vector s1 (no preconditioning)
! -----------------------------------------------------------------

module test_fsunlinsol_spgmr_serial
  use, intrinsic :: iso_c_binding
  use test_utilities
  implicit none

  integer(kind=myindextype), private, parameter :: N = 100
  integer(c_int), private, parameter :: pretype = 1     ! Preconditioning type (1 or 2)
  integer(c_int), private, parameter :: gstype = 1     ! Gram-Schmidt orthoognalization type (1 or 2)
  integer(c_int), private, parameter :: maxl = 100   ! maximum Krylov subspace dimension (> 0)
  real(c_double), private, parameter :: tol = 1e-13 ! solver tolerance

  type, private :: UserData
    integer(kind=myindextype) :: N
    type(N_Vector), pointer  :: d, s1, s2
  end type

contains

  integer(c_int) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod
    use fsunlinsol_spgmr_mod
    use test_sunlinsol

    implicit none

    type(SUNLinearSolver), pointer :: LS         ! test linear solver
    type(SUNMatrix), pointer :: A          ! dummy SUNMatrix
    type(N_Vector), pointer :: x, xhat, b ! test vectors
    type(UserData), pointer :: probdata   ! problem data
    real(c_double), pointer :: xdata(:)   ! x vector data
    real(c_double)                 :: tmpr       ! temporary real value
    integer(kind=myindextype)     :: j
    integer(c_int)                 :: tmp

    ! setup
    fails = 0

    A => null()

    x => FN_VNew_Serial(N, sunctx)
    xhat => FN_VNew_Serial(N, sunctx)
    b => FN_VNew_Serial(N, sunctx)

    allocate (probdata)
    probdata%N = N
    probdata%d => FN_VNew_Serial(N, sunctx)
    probdata%s1 => FN_VNew_Serial(N, sunctx)
    probdata%s2 => FN_VNew_Serial(N, sunctx)

    ! fill xhat vector with uniform random data in [1, 2)
    xdata => FN_VGetArrayPointer(xhat)
    do j = 1, N
      call random_number(tmpr)
      xdata(j) = ONE + tmpr
    end do

    ! fill Jacobi vector with matrix diagonal
    call FN_VConst(FIVE, probdata%d)

    ! create SPGMR linear solver
    LS => FSUNLinSol_SPGMR(x, pretype, maxl, sunctx)

    ! run initialization tests
    fails = fails + Test_FSUNLinSolGetType(LS, SUNLINEARSOLVER_ITERATIVE, 0)
    fails = fails + Test_FSUNLinSolSetATimes(LS, c_loc(probdata), &
                                             c_funloc(ATimes), 0)
    fails = fails + Test_FSUNLinSolSetPreconditioner(LS, &
                                                     c_loc(probdata), &
                                                     c_funloc(PSetup), &
                                                     c_funloc(PSolve), &
                                                     0)
    fails = fails + Test_FSUNLinSolSetScalingVectors(LS, probdata%s1, &
                                                     probdata%s2, 0)
    fails = fails + Test_FSUNLinSolInitialize(LS, 0)
    fails = fails + Test_FSUNLinSolSpace(LS, 0)
    fails = fails + FSUNLinSol_SPGMRSetGSType(LS, gstype)

    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_SPGMR module, initialization'
    else
      print *, 'SUCCESS: FSUNLinSol_SPGMR module, initialization tests passed'
      print *, ''
    end if

    ! Test 1: simple Poisson-like solve (no preconditioning)

    ! set scaling vectors
    call FN_VConst(ONE, probdata%s1)
    call FN_VConst(ONE, probdata%s2)

    ! fill x vector with scaled version
    call FN_VDiv(xhat, probdata%s2, x)

    ! fill b vector with result of matrix-vector product
    fails = fails + ATimes(c_loc(probdata), x, b)

    ! Run tests with this setup
    fails = fails + FSUNLinSol_SPGMRSetPrecType(LS, SUN_PREC_NONE); 
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0); 
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, tol, 0); 
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0); 
    fails = fails + Test_FSUNLinSolNumIters(LS, 0); 
    fails = fails + Test_FSUNLinSolResNorm(LS, 0); 
    fails = fails + Test_FSUNLinSolResid(LS, 0); 
    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_SPGMR module, problem 1'
    else
      print *, 'SUCCESS: FSUNLinSol_SPGMR module, problem 1, all tests passed'
      print *, ''
    end if

    ! Test 2: simple Poisson-like solve (Jacobi preconditioning)

    ! set scaling vectors
    call FN_VConst(ONE, probdata%s1)
    call FN_VConst(ONE, probdata%s2)

    ! fill x vector with scaled version
    call FN_VDiv(xhat, probdata%s2, x)

    ! fill b vector with result of matrix-vector product
    fails = fails + ATimes(c_loc(probdata), x, b)

    ! Run tests with this setup
    fails = fails + FSUNLinSol_SPGMRSetPrecType(LS, pretype); 
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0); 
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, tol, 0); 
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0); 
    fails = fails + Test_FSUNLinSolNumIters(LS, 0); 
    fails = fails + Test_FSUNLinSolResNorm(LS, 0); 
    fails = fails + Test_FSUNLinSolResid(LS, 0); 
    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_SPGMR module, problem 2'
    else
      print *, 'SUCCESS: FSUNLinSol_SPGMR module, problem 2, all tests passed'
      print *, ''
    end if

    ! Test 3: Poisson-like solve w/ scaled rows (no preconditioning)

    ! set scaling vectors
    xdata => FN_VGetArrayPointer(probdata%s1)
    do j = 1, N
      call random_number(tmpr)
      xdata(j) = ONE + 1000.0d0*tmpr
    end do
    call FN_VConst(ONE, probdata%s2)

    ! fill x vector with scaled version
    call FN_VDiv(xhat, probdata%s2, x)

    ! fill b vector with result of matrix-vector product
    fails = fails + ATimes(c_loc(probdata), x, b)

    ! Run tests with this setup
    fails = fails + FSUNLinSol_SPGMRSetPrecType(LS, SUN_PREC_NONE); 
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0); 
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, tol, 0); 
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0); 
    fails = fails + Test_FSUNLinSolNumIters(LS, 0); 
    fails = fails + Test_FSUNLinSolResNorm(LS, 0); 
    fails = fails + Test_FSUNLinSolResid(LS, 0); 
    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_SPGMR module, problem 3'
    else
      print *, 'SUCCESS: FSUNLinSol_SPGMR module, problem 3, all tests passed'
      print *, ''
    end if

    ! cleanup
    tmp = FSUNLinSolFree(LS)
    call FN_VDestroy(x)
    call FN_VDestroy(xhat)
    call FN_VDestroy(b)
    call FN_VDestroy(probdata%d)
    call FN_VDestroy(probdata%s1)
    call FN_VDestroy(probdata%s2)
    deallocate (probdata)

  end function unit_tests

  integer(c_int) function ATimes(udata, vvec, zvec) result(ret) bind(C)
    use, intrinsic :: iso_c_binding

    use test_utilities

    implicit none

    type(c_ptr), value :: udata
    type(N_Vector) :: vvec, zvec
    type(UserData), pointer :: probdata
    real(c_double), pointer :: v(:), z(:), s1(:), s2(:)
    integer(c_long) :: i, N

    call c_f_pointer(udata, probdata)

    v => FN_VGetArrayPointer(vvec)
    z => FN_VGetArrayPointer(zvec)
    s1 => FN_VGetArrayPointer(probdata%s1)
    s2 => FN_VGetArrayPointer(probdata%s2)
    N = probdata%N

    ! perform product at the left domain boundary (note: v is zero at the boundary)
    z(1) = (FIVE*v(1)*s2(1) - v(2)*s2(2))/s1(1)

    ! iterate through interior of local domain, performing product
    do i = 2, N - 1
      z(i) = (-v(i - 1)*s2(i - 1) + FIVE*v(i)*s2(i) - v(i + 1)*s2(i + 1))/s1(i)
    end do

    ! perform product at the right domain boundary (note: v is zero at the boundary)
    z(N) = (-v(N - 1)*s2(N - 1) + FIVE*v(N)*s2(N))/s1(N)

    ret = 0
  end function ATimes

  integer(c_int) function PSetup(udata) result(ret) bind(C)
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: udata
    ret = 0
  end function PSetup

  integer(c_int) function PSolve(udata, rvec, zvec, tol, lr) &
    result(ret) bind(C)
    use, intrinsic :: iso_c_binding

    use test_utilities

    implicit none

    type(c_ptr), value      :: udata
    type(N_Vector)             :: rvec, zvec
    real(c_double)             :: tol
    integer(c_int)             :: lr
    type(UserData), pointer    :: probdata
    real(c_double), pointer    :: r(:), z(:), d(:)
    integer(kind=myindextype) :: i, N

    call c_f_pointer(udata, probdata)

    r => FN_VGetArrayPointer(rvec)
    z => FN_VGetArrayPointer(zvec)
    d => FN_VGetArrayPointer(probdata%d)
    N = probdata%N

    do i = 1, N
      z(i) = r(i)/d(i)
    end do

    ret = 0
  end function PSolve

end module

integer(c_int) function check_vector(X, Y, tol) result(failure)
  use, intrinsic :: iso_c_binding
  use test_fsunlinsol_spgmr_serial
  use test_utilities
  implicit none

  type(N_Vector)  :: x, y
  real(c_double)  :: tol, maxerr
  integer(kind=myindextype) :: i, xlen, ylen
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
    failure = failure + FNEQTOL(xdata(i), ydata(i), FIVE*tol*abs(xdata(i)))
  end do

  if (failure > 0) then
    maxerr = ZERO
    do i = 1, xlen
      maxerr = max(abs(xdata(i) - ydata(i))/abs(xdata(i)), maxerr)
    end do
    write (*, '(A,E14.7,A,E14.7,A)') &
      "FAIL: check_vector failure: maxerr = ", maxerr, "  (tol = ", FIVE*tol, ")"
  end if

end function check_vector

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunlinsol_spgmr_serial

  !======== Declarations ========
  implicit none
  integer(c_int) :: fails = 0

  !============== Introduction =============
  print *, 'SPGMR SUNLinearSolver Fortran 2003 interface test'
  print *, ''

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
