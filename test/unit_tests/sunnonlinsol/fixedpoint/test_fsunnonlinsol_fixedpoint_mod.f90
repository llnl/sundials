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
! fixedpoint SUNNonlinearSolver implementation.
! -----------------------------------------------------------------

module test_fsunnonlinsol_fixedpoint
  use, intrinsic :: iso_c_binding
  use test_utilities

  implicit none

  integer(kind=myindextype), parameter :: NEQ = 3 ! number of equations
  integer(c_int), parameter :: MAXIT = 20     ! max nonlinear iters.
  real(c_double), parameter :: TOL = 1.0e-4 ! nonlinear solver tolerance

  real(c_double), parameter :: PI = 3.1415926535898

  ! approximate solution
  real(c_double) :: XTRUE = 0.5d0
  real(c_double) :: YTRUE = 1.0d0
  real(c_double) :: ZTRUE = -PI/6.0d0

  type(N_Vector), pointer :: y0

contains

  integer(c_int) function unit_tests() result(retval)
    use, intrinsic :: iso_c_binding
    use fsundials_core_mod
    use fnvector_serial_mod
    use fsunnonlinsol_fixedpoint_mod

    implicit none

    type(SUNNonlinearSolver), pointer :: NLS         ! test nonlinear solver
    type(N_Vector), pointer :: ycur, ycor, w ! test vectors
    real(c_double), pointer :: data(:)
    integer(c_long)                   :: niters(1)
    integer(c_int)                    :: tmp

    y0 => FN_VNew_Serial(NEQ, sunctx)
    ycor => FN_VClone(y0)
    ycur => FN_VClone(y0)
    w => FN_VClone(y0)

    ! set initial guess
    data => FN_VGetArrayPointer(y0)
    data(1) = 0.1d0
    data(2) = 0.1d0
    data(3) = -0.1d0

    ! set initial correction
    call FN_VConst(0.0d0, ycor)

    ! set weights
    call FN_VConst(1.0d0, w)

    ! create and test NLS
    NLS => FSUNNonlinsol_FixedPoint(y0, 0, sunctx)

    retval = FSUNNonlinSolSetSysFn(NLS, c_funloc(FPFunction))
    if (retval /= 0) then
      write (*, '(A,I0)') '   >>> FAIL: FSUNNonlinSolSetSysFn returned ', retval
      return
    end if

    retval = FSUNNonlinSolSetConvTestFn(NLS, c_funloc(ConvTest), c_null_ptr)
    if (retval /= 0) then
      write (*, '(A,I0)') '   >>> FAIL: FSUNNonlinSolSetConvTestFn returned ', retval
      return
    end if

    retval = FSUNNonlinSolSetMaxIters(NLS, MAXIT)
    if (retval /= 0) then
      write (*, '(A,I0)') '   >>> FAIL: FSUNNonlinSolSetMaxIters returned ', retval
      return
    end if

    retval = FSUNNonlinSolSolve(NLS, y0, ycor, w, TOL, 1, c_loc(y0))
    if (retval /= 0) then
      write (*, '(A,I0)') '   >>> FAIL: FSUNNonlinSolSolve returned ', retval
      return
    end if

    ! update the initial guess with the final correction
    call FN_VLinearSum(1.0d0, y0, 1.0d0, ycor, ycur); 
    ! print number of iterations
    retval = FSUNNonlinSolGetNumIters(NLS, niters)
    if (retval /= 0) then
      write (*, '(A,I0)') '   >>> FAIL: FSUNNonlinSolGetNumIters returned ', retval
      return
    end if

    write (*, '(A,I0)') 'Number of nonlinear iterations: ', niters(1)

    ! check answer
    retval = check_ans(ycur, TOL)
    if (retval /= 0) then
      write (*, '(A,I0)') '   >>> FAIL: check_ans failed'
      return
    end if

    ! cleanup
    call FN_VDestroy(y0)
    call FN_VDestroy(ycor)
    call FN_VDestroy(ycur)
    call FN_VDestroy(w)
    tmp = FSUNNonlinSolFree(NLS)

  end function unit_tests

  integer(c_int) function ConvTest(NLS, y, del, tol, ewt, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNNonlinearSolver) :: NLS
    type(N_Vector)           :: y, del, ewt
    real(c_double), value    :: tol
    type(c_ptr), value       :: mem
    real(c_double)           :: delnrm

    ! compute the norm of the correction
    delnrm = FN_VMaxNorm(del)

    if (delnrm <= tol) then
      retval = SUN_SUCCESS  ! converged
    else
      retval = SUN_NLS_CONTINUE ! not converged
    end if

  end function

  ! -----------------------------------------------------------------------------
  ! Nonlinear system F(x,y,z):
  !
  ! 3x - cos((y-1)z) - 1/2 = 0
  ! x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
  ! exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
  !
  ! Nonlinear fixed point function G(x,y,z):
  !
  ! G1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
  ! G2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
  ! G3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
  !
  ! Corrector form g(x,y,z):
  !
  ! g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6 - x0
  ! g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9 - y0
  ! g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60 - z0
  !
  ! ----------------------------------------------------------------------------
  integer(c_int) function FPFunction(ycor, f, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    implicit none

    type(N_Vector)          :: ycor, f
    type(c_ptr), value      :: mem
    real(c_double), pointer :: data(:), fdata(:)
    real(c_double)          :: x, y, z

    data => FN_VGetArrayPointer(ycor)
    fdata => FN_VGetArrayPointer(f)

    x = data(1)
    y = data(2)
    z = data(3)

    fdata(1) = (1.0d0/3.0d0)*cos((y - 1.0d0)*z) + (1.0d0/6.0d0)
    fdata(2) = (1.0d0/9.0d0)*sqrt(x*x + sin(z) + 1.06d0) + 0.9d0
    fdata(3) = -(1/20.d0)*exp(-x*(y - 1.0d0)) - (10.d0*PI - 3.0d0)/60.0d0

    call FN_VLinearSum(1.0d0, f, -1.0d0, y0, f)

    retval = 0

  end function

  integer(c_int) function check_ans(ycor, tol) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    implicit none

    type(N_Vector) :: ycor
    real(c_double), value :: tol
    real(c_double) ::  ex, ey, ez
    real(c_double), pointer :: data(:)

    ! extract and print solution
    data => FN_VGetArrayPointer(ycor)

    write (*, *) 'Solution:'
    write (*, '(A,E14.7)') '    x = ', data(1)
    write (*, '(A,E14.7)') '    y = ', data(2)
    write (*, '(A,E14.7)') '    z = ', data(3)

    ex = data(1) - XTRUE
    ey = data(2) - YTRUE
    ez = data(3) - ZTRUE

    write (*, *) 'Solution Error:'
    write (*, '(A,E14.7)') '    ex = ', ex
    write (*, '(A,E14.7)') '    ey = ', ey
    write (*, '(A,E14.7)') '    ez = ', ez

    tol = tol*10.0d0
    if (ex > tol .or. ey > tol .or. ez > tol) then
      retval = 1
    end if

    retval = 0
  end function

end module

program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunnonlinsol_fixedpoint

  !======== Declarations ========
  implicit none

  integer(c_int) :: fails = 0

  !============== Introduction =============
  write (*, *) 'SUNNonlinearSolver_FixedPoint Fortran 2003 interface test'

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
