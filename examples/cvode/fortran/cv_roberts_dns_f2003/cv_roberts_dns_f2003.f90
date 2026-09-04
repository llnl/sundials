! ------------------------------------------------------------------
! Programmer(s): Daniel M. Margolis @ UMBC
!                based off the previous Fortran-77 example program,
!                cvode/fcmix_serial/fcvRoberts_dns.f
! ------------------------------------------------------------------
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
! ------------------------------------------------------------------
! The following is a simple example problem for CVODE, due to Robertson,
! is from chemical kinetics, and consists of the following three rate
! equations:
!
!      dy1/dt = -.04*y1 + 1.e4*y2*y3
!      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!      dy3/dt = 3.e7*y2**2
!
! on the interval from t = 0.0 to t = 4.e10, with initial
! conditions: y1 = 1, y2 = y3 = 0.
!
! While integrating the system, we also use the rootfinding
! feature to find the points at which y1 = 1.e-4 or at which
! y3 = 0.01.
!
! The problem is solved with CVODE using the DENSE linear
! solver, with a user-supplied Jacobian. Output is printed at
! t = .4, 4, 40, ..., 4e10. It uses ATOL much smaller for y2
! than y1 or y3 because y2 has much smaller values. At the end
! of the run, various counters of interest are printed.
! ------------------------------------------------------------------

module robertsDns_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  !======= Declarations =========
  implicit none

  integer(c_int), parameter :: nout = 12
  integer(c_int64_t), parameter :: neq = 3

contains

  ! ----------------------------------------------------------------
  ! fcnrob: The CVODE RHS operator function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function fcnrob(t, sunvec_y, sunvec_f, user_data) &
    result(ierr) bind(C, name='fcnrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! function N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yval(:)
    real(c_double), pointer, dimension(neq) :: fval(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    fval => FN_VGetArrayPointer(sunvec_f)

    ! fill residual vector
    fval(1) = -0.04d0*yval(1) + 1.0d4*yval(2)*yval(3)
    fval(3) = 3.0d7*yval(2)**2
    fval(2) = -fval(1) - fval(3)

    ! return success
    ierr = 0
    return

  end function fcnrob
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! grob: The root function routine
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function grob(t, sunvec_y, gout, user_data) &
    result(ierr) bind(C, name='grob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    real(c_double)        :: gout(2)   ! root function values
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yval(:)

    !======= Internals ============

    ! get data array from SUNDIALS vector
    yval => FN_VGetArrayPointer(sunvec_y)

    ! fill root vector
    gout(1) = yval(1) - 1.0d-4
    gout(2) = yval(3) - 1.0d-2

    ! return success
    ierr = 0
    return

  end function grob
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! jacrob: The Jacobian function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function jacrob(t, sunvec_y, sunvec_f, &
                                 sunmat_J, user_data, sunvec_t1, sunvec_t2, sunvec_t3) &
    result(ierr) bind(C, name='jacrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsunmatrix_dense_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! residual N_Vector
    type(SUNMatrix)       :: sunmat_J  ! Jacobian SUNMatrix
    type(c_ptr), value :: user_data ! user-defined data
    type(N_Vector)        :: sunvec_t1 ! temporary N_Vectors
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! pointers to data in SUNDIALS vector and matrix
    real(c_double), pointer, dimension(neq) :: yval(:)
    real(c_double), pointer, dimension(neq, neq) :: J(:, :)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    J(1:3, 1:3) => FSUNDenseMatrix_Data(sunmat_J)

    ! fill Jacobian entries
    J(1, 1) = -0.04d0
    J(2, 1) = 0.04d0
    J(3, 1) = 0.0d0
    J(1, 2) = 1.0d4*yval(3)
    J(2, 2) = -1.0d4*yval(3) - 6.0d7*yval(2)
    J(3, 2) = 6.0d7*yval(2)
    J(1, 3) = 1.0d4*yval(2)
    J(2, 3) = -1.0d4*yval(2)
    J(3, 3) = 0.0d0

    ! return success
    ierr = 0
    return

  end function jacrob
  ! ----------------------------------------------------------------

end module robertsDns_mod
! ------------------------------------------------------------------

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fcvode_mod                    ! Fortran interface to CVODE
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
  use robertsDns_mod                ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: rtol, t0, tout, tret(1)
  integer(c_int) :: iout, retval, retvalr, nrtfn, rootsfound(2)

  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_dky    ! sundials solution vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(c_ptr)                       :: cvode_mem     ! CVode memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context

  ! solution and tolerance vectors, neq is set in the robertsDns_mod module
  real(c_double) :: yval(neq), avtol(neq), dkyval(neq)

  !======= Internals ============

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! initialize solution vectors and tolerances
  yval(1) = 1.0d0
  yval(2) = 0.0d0
  yval(3) = 0.0d0

  rtol = 1.0d-4

  avtol(1) = 1.0d-8
  avtol(2) = 1.0d-14
  avtol(3) = 1.0d-6

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(neq, yval, sunctx)
  if (.not. associated(sunvec_y)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  sunvec_av => FN_VMake_Serial(neq, avtol, sunctx)
  if (.not. associated(sunvec_av)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  ! set limits
  t0 = 0.0d0
  tout = 0.4d0

  call PrintHeader(rtol, avtol, yval)

  ! Call FCVodeCreate and FCVodeInit to create and initialize CVode memory
  cvode_mem = FCVodeCreate(CV_BDF, sunctx)
  if (.not. c_associated(cvode_mem)) print *, 'ERROR: cvode_mem = NULL'

  retval = FCVodeInit(cvode_mem, c_funloc(fcnrob), t0, sunvec_y)
  if (retval /= 0) then
    print *, 'Error in FCVodeInit, retval = ', retval, '; halting'
    stop 1
  end if

  ! Call FCVodeSVtolerances to set tolerances
  retval = FCVodeSVtolerances(cvode_mem, rtol, sunvec_av)
  if (retval /= 0) then
    print *, 'Error in FCVodeSVtolerances, retval = ', retval, '; halting'
    stop 1
  end if

  ! Call FCVodeRootInit to specify the root function grob with 2 components
  nrtfn = 2
  retval = FCVodeRootInit(cvode_mem, nrtfn, c_funloc(grob))
  if (retval /= 0) then
    print *, 'Error in FCVodeRootInit, retval = ', retval, '; halting'
    stop 1
  end if

  ! Create dense SUNMatrix for use in linear solves
  sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
  if (.not. associated(sunmat_A)) then
    print *, 'ERROR: sunmat = NULL'
    stop 1
  end if

  ! Create dense SUNLinearSolver object
  sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
  if (.not. associated(sunlinsol_LS)) then
    print *, 'ERROR: sunlinsol = NULL'
    stop 1
  end if

  ! Attach the matrix and linear solver
  retval = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A); 
  if (retval /= 0) then
    print *, 'Error in FCVodeSetLinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! Set the user-supplied Jacobian routine
  retval = FCVodeSetJacFn(cvode_mem, c_funloc(jacrob))
  if (retval /= 0) then
    print *, 'Error in FCVodeSetJacFn, retval = ', retval, '; halting'
    stop 1
  end if

  ! In loop, call FCVode, print results, and test for error.

  iout = 0
  do while (iout < nout)

    retval = FCVode(cvode_mem, tout, sunvec_y, tret(1), CV_NORMAL)
    if (retval < 0) then
      print *, 'Error in FCVode, retval = ', retval, '; halting'
      stop 1
    end if

    call PrintOutput(cvode_mem, tret(1), yval)

    if (retval == CV_ROOT_RETURN) then
      retvalr = FCVodeGetRootInfo(cvode_mem, rootsfound)
      if (retvalr < 0) then
        print *, 'Error in FCVodeGetRootInfo, retval = ', retval, '; halting'
        stop 1
      end if
      print '(a,2(i2,2x))', "    rootsfound[] = ", rootsfound(1), rootsfound(2)
    end if

    if (retval == CV_SUCCESS) then
      iout = iout + 1
      tout = tout*10.0d0
    end if
  end do

  sunvec_dky => FN_VMake_Serial(neq, dkyval, sunctx)
  if (.not. associated(sunvec_dky)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  ! find and print derivative at tret(1)
  retval = FCVodeGetDky(cvode_mem, tret(1), 1, sunvec_dky)
  if (retval /= 0) then
    print *, 'Error in CVodeGetDky'
    stop 1
  end if
  print *, " "
  print *, "---------------------------------------------------"
  print *, "     Final      y1'          y2'          y3'"
  print *, "---------------------------------------------------"
  print '(13x,3(es12.4,1x))', dkyval

  call PrintFinalStats(cvode_mem)

  ! free memory
  call FCVodeFree(cvode_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_dky)
  call FN_VDestroy(sunvec_av)
  retval = FSUNContext_Free(sunctx)

end program main
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintHeader: prints first lines of output (problem description)
! ----------------------------------------------------------------
subroutine PrintHeader(rtol, avtol, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use robertsDns_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double) :: rtol
  real(c_double) :: avtol(neq)
  real(c_double) :: y(neq)

  !======= Internals ============

  print *, " "
  print *, "cv_roberts_dns_f2003.f90: Robertson CV ODE serial example problem for CVODE"
  print *, "         Three equation chemical kinetics problem."
  print *, " "
  print *, "Linear solver: DENSE, with user-supplied Jacobian."
  print '(a,f6.4,a,3(es7.0,1x))', "Tolerance parameters:  rtol = ", rtol, "   atol = ", avtol
  print '(a,3(f5.2,1x),a)', "Initial conditions y0 = (", y, ")"
  print *, "Constraints not used."
  print *, " "
  print *, "---------------------------------------------------"
  print *, "   t            y1           y2           y3"
  print *, "---------------------------------------------------"

  return
end subroutine PrintHeader
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintOutput
! ----------------------------------------------------------------
subroutine PrintOutput(cvode_mem, t, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use robertsDns_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  type(c_ptr)    :: cvode_mem
  real(c_double) :: t, y(neq)

  !======= Internals ============

  print '(es12.4,1x,3(es12.4,1x))', t, y(1), y(2), y(3)

end subroutine PrintOutput
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintFinalStats
!
! Print ARKSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(cvode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use fcvode_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

  integer(c_int)  :: retval          ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nfe(1)        ! num function evals
  integer(c_long) :: netfails(1)   ! num error test fails
  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nncfails(1)   ! nonlinear solver fails
  integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations
  integer(c_long) :: nluevals(1)   ! number of LU evals
  integer(c_long) :: ngevals(1)    ! number of root evals

  !======= Internals ============

  retval = FCVodeGetNumSteps(cvode_mem, nsteps)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumSteps, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumRhsEvals(cvode_mem, nfe)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumRhsEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumLinSolvSetups(cvode_mem, nluevals)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumLinSolvSetups, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumErrTestFails(cvode_mem, netfails)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumErrTestFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumNonlinSolvIters(cvode_mem, nniters)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumNonlinSolvIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumNonlinSolvConvFails(cvode_mem, nncfails)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumJacEvals(cvode_mem, njacevals)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumJacEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeGetNumGEvals(cvode_mem, ngevals)
  if (retval /= 0) then
    print *, 'Error in FCVodeGetNumGEvals, retval = ', retval, '; halting'
    stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)', 'Total internal steps taken    =', nsteps
  print '(4x,A,i9)', 'Total rhs function calls      =', nfe
  print '(4x,A,i9)', 'Total Jacobian function calls =', njacevals
  print '(4x,A,i9)', 'Total root function calls     =', ngevals
  print '(4x,A,i9)', 'Total LU function calls       =', nluevals
  print '(4x,A,i9)', 'Num error test failures       =', netfails
  print '(4x,A,i9)', 'Num nonlinear solver iters    =', nniters
  print '(4x,A,i9)', 'Num nonlinear solver fails    =', nncfails
  print *, ' '

  return

end subroutine PrintFinalStats
! ----------------------------------------------------------------
