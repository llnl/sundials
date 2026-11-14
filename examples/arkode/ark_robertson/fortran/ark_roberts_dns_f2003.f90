! ------------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ UMBC
!                modified by Daniel M. Margolis @ UMBC
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
! The following is a simple example problem for ARKODE, due to Robertson,
! is from chemical kinetics, and consists of the following three
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
! feature to find the points at which y1 = 1e-4 or at which
! y3 = 0.01.
!
! The problem is solved with ARKODE using the DENSE linear
! solver, with a user-supplied Jacobian. Output is printed at
! t = .4, 4, 40, ..., 4e10.
! ------------------------------------------------------------------

module dns_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  !======= Declarations =========
  implicit none

  integer(c_int64_t), parameter :: neq = 3
  integer(c_long), parameter :: nout = 12

contains

  ! ----------------------------------------------------------------
  ! fcnirob: The implicit RHS operator function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function fcnirob(tn, sunvec_y, sunvec_f, user_data) &
    result(ierr) bind(C, name='fcnirob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
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

  end function fcnirob
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! grob: The root function routine
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function grob(tn, sunvec_y, gout, user_data) &
    result(ierr) bind(C, name='grob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    real(c_double)        :: gout(2)   ! root function values
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yval(:)

    !======= Internals ============

    ! get data array from SUNDIALS vector
    yval => FN_VGetArrayPointer(sunvec_y)

    ! fill root vector
    gout(1) = yval(1) - 0.0001d0
    gout(2) = yval(3) - 0.01d0

    ! return success
    ierr = 0
    return

  end function grob
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! jacrob: The ODE Jacobian function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function jacrob(tn, sunvec_y, sunvec_f, &
                                 sunmat_J, user_data, sunvec_t1, sunvec_t2, sunvec_t3) &
    result(ierr) bind(C, name='jacrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod
    use fsunmatrix_dense_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! residual N_Vector
    type(SUNMatrix)       :: sunmat_J  ! Jacobian SUNMatrix
    type(c_ptr), value :: user_data ! user-defined data
    type(N_Vector)        :: sunvec_t1 ! temporary N_Vectors
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! pointers to data in SUNDIALS vector and matrix
    real(c_double), pointer, dimension(neq)     :: yval(:)
    real(c_double), pointer, dimension(neq, neq) :: J(:, :)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    J(1:3, 1:3) => FSUNDenseMatrix_Data(sunmat_J)

    ! fill Jacobian entries
    J(1, 1) = -0.04d0
    J(2, 1) = 0.04d0
    J(3, 1) = 0.d0
    J(1, 2) = 1.d4*yval(3)
    J(2, 2) = -1.d4*yval(3) - 6.0d7*yval(2)
    J(3, 2) = 6.d7*yval(2)
    J(1, 3) = 1.d4*yval(2)
    J(2, 3) = -1.d4*yval(2)
    J(3, 3) = 0.d0

    ! return success
    ierr = 0
    return

  end function jacrob
  ! ----------------------------------------------------------------

end module dns_mod
! ------------------------------------------------------------------

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                   ! Fortran interface to ARKODE
  use farkode_arkstep_mod           ! Fortran interface to the ARKStep module
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
  use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver
  use dns_mod                       ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: rtol, t0, tout1, tout, tret(1)
  integer(c_int) :: iout, retval, retvalr, nrtfn, rootsfound(2)

  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_dky    ! sundials solution vector
  type(N_Vector), pointer :: sunvec_f      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: arkode_mem    ! ARKODE memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context

  ! solution and tolerance vectors, neq is set in the dns_mod module
  real(c_double) :: yval(neq), fval(neq), avtol(neq), dkyval(neq)

  ! fine-tuning initialized here
  real(c_double)  :: initsize, nlscoef
  integer(c_long) :: mxsteps
  integer(c_int)  :: nliters, pmethod, maxetf

  !======= Internals ============

  ! create the SUNDIALS context
  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! initialize solution vectors and tolerances
  yval(1) = 1.d0
  yval(2) = 0.d0
  yval(3) = 0.d0
  fval = 0.d0
  rtol = 1.d-4
  avtol(1) = 1.d-8
  avtol(2) = 1.d-11
  avtol(3) = 1.d-8

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(neq, yval, sunctx)
  if (.not. associated(sunvec_y)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  sunvec_f => FN_VMake_Serial(neq, fval, sunctx)
  if (.not. associated(sunvec_f)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  sunvec_av => FN_VMake_Serial(neq, avtol, sunctx)
  if (.not. associated(sunvec_av)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  ! set integration limits
  t0 = 0.d0
  tout1 = 0.4d0

  call PrintHeader(rtol, avtol, yval)

  ! Call FARKStepCreate to initialize ARKODE memory
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(fcnirob), t0, sunvec_y, sunctx)
  if (.not. c_associated(arkode_mem)) print *, 'ERROR: arkode_mem = NULL'

  ! Call FARKodeSVtolerances to set tolerances
  retval = FARKodeSVtolerances(arkode_mem, rtol, sunvec_av)
  if (retval /= 0) then
    print *, 'Error in FARKodeSVtolerances, retval = ', retval, '; halting'
    stop 1
  end if

  ! Call FARKodeRootInit to specify the root function grob with 2 components
  nrtfn = 2
  retval = FARKodeRootInit(arkode_mem, nrtfn, c_funloc(grob))
  if (retval /= 0) then
    print *, 'Error in FARKodeRootInit, retval = ', retval, '; halting'
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
  retval = FARKodeSetLinearSolver(arkode_mem, sunlinsol_LS, sunmat_A); 
  if (retval /= 0) then
    print *, 'Error in FARKodeSetLinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! Set the user-supplied Jacobian routine
  retval = FARKodeSetJacFn(arkode_mem, c_funloc(jacrob))
  if (retval /= 0) then
    print *, 'Error in FARKodeSetJacFn, retval = ', retval, '; halting'
    stop 1
  end if

  ! Set additional method parameters
  mxsteps = 10000
  retval = FARKodeSetMaxNumSteps(arkode_mem, mxsteps)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetMaxNumSteps'
    stop 1
  end if

  initsize = 1.d-4*rtol
  retval = FARKodeSetInitStep(arkode_mem, initsize)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetInitStep'
    stop 1
  end if

  nlscoef = 1.d-7
  retval = FARKodeSetNonlinConvCoef(arkode_mem, nlscoef)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetNonlinConvCoef'
    stop 1
  end if

  nliters = 8
  retval = FARKodeSetMaxNonlinIters(arkode_mem, nliters)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetMaxNonlinIters'
    stop 1
  end if

  pmethod = 1
  retval = FARKodeSetPredictorMethod(arkode_mem, pmethod)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetPredictorMethod'
    stop 1
  end if

  maxetf = 20
  retval = FARKodeSetMaxErrTestFails(arkode_mem, maxetf)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetMaxErrTestFails'
    stop 1
  end if

  ! Create Newton SUNNonlinearSolver object. ARKODE uses a
  ! Newton SUNNonlinearSolver by default, so it is not necessary
  ! to create it and attach it. It is done in this example code
  ! solely for demonstration purposes.
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
  if (.not. associated(sunnonlin_NLS)) then
    print *, 'ERROR: sunnonlinsol = NULL'
    stop 1
  end if

  ! Attach the nonlinear solver
  retval = FARKodeSetNonlinearSolver(arkode_mem, sunnonlin_NLS)
  if (retval /= 0) then
    print *, 'Error in FARKodeSetNonlinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! In loop, call ARKodeEvolve, print results, and test for error.
  iout = 0
  tout = tout1
  do while (iout < nout)

    retval = FARKodeEvolve(arkode_mem, tout, sunvec_y, tret(1), ARK_NORMAL)
    if (retval < 0) then
      print *, 'Error in FARKodeEvolve, retval = ', retval, '; halting'
      stop 1
    end if

    call PrintOutput(arkode_mem, tret(1), yval)

    if (retval == ARK_ROOT_RETURN) then
      retvalr = FARKodeGetRootInfo(arkode_mem, rootsfound)
      if (retvalr < 0) then
        print *, 'Error in FARKodeGetRootInfo, retval = ', retval, '; halting'
        stop 1
      end if
      print '(a,2(i2,2x))', "    rootsfound[] = ", rootsfound(1), rootsfound(2)
    end if

    if (retval == ARK_SUCCESS) then
      iout = iout + 1
      tout = tout*10.d0
    end if
  end do

  ! find and print derivative at tret(1)
  sunvec_dky => FN_VMake_Serial(neq, dkyval, sunctx)
  if (.not. associated(sunvec_dky)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  retval = FARKodeGetDky(arkode_mem, tret(1), 1, sunvec_dky)
  if (retval /= 0) then
    print *, 'Error in ARKodeGetDky'
    stop 1
  end if
  print *, " "
  print *, "------------------------------------------------------"
  print *, "     Final      y1'          y2'          y3'"
  print *, "------------------------------------------------------"
  print '(13x,3(es12.4,1x))', dkyval

  ! diagnostics output
  call PrintFinalStats(arkode_mem)

  ! free memory
  call FARKodeFree(arkode_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_f)
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
  use dns_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double) :: rtol
  real(c_double) :: avtol(neq)
  real(c_double) :: y(neq)

  !======= Internals ============

  print *, " "
  print *, "ark_roberts_dns_f2003.f90: Robertson ARK ODE serial example problem for ARKODE"
  print *, "         Three equation chemical kinetics problem."
  print *, " "
  print *, "Linear solver: DENSE, with user-supplied Jacobian."
  print '(a,f6.4,a,3(es7.0,1x))', "Tolerance parameters:  rtol = ", rtol, "   atol = ", avtol
  print '(a,3(f5.2,1x),a)', "Initial conditions y0 = (", y, ")"
  print *, "Constraints not used."
  print *, " "
  print *, "----------------------------------------------------------------------"
  print *, "   t            y1           y2           y3       | nst      h"
  print *, "----------------------------------------------------------------------"

  return
end subroutine PrintHeader
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintOutput
! ----------------------------------------------------------------
subroutine PrintOutput(arkode_mem, t, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod
  use dns_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  type(c_ptr)    :: arkode_mem
  real(c_double) :: t, y(neq)

  ! internal variables
  integer(c_int)  :: retval
  integer(c_long) :: nst(1)
  real(c_double)  :: hused(1)

  !======= Internals ============

  retval = FARKodeGetNumSteps(arkode_mem, nst)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumSteps, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetLastStep(arkode_mem, hused)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetLastStep, retval = ', retval, '; halting'
    stop 1
  end if

  print '(es12.4,1x,3(es12.4,1x),a,i3,2x,es12.4)', &
    t, y(1), y(2), y(3), "| ", nst, hused(1)

end subroutine PrintOutput
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintFinalStats
!
! Print ARKStep statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(arkode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: arkode_mem ! solver memory structure

  integer(c_int)  :: retval          ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nst_a(1)      ! num steps attempted
  integer(c_long) :: nfe(1)        ! num explicit function evals
  integer(c_long) :: nfi(1)        ! num implicit function evals
  integer(c_long) :: nlinsetups(1) ! num linear solver setups
  integer(c_long) :: netfails(1)   ! num error test fails

  real(c_double)  :: hinused(1)    ! initial step size
  real(c_double)  :: hlast(1)      ! last step size
  real(c_double)  :: hcur(1)       ! step size for next step
  real(c_double)  :: tcur(1)       ! internal time reached

  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nncfails(1)   ! nonlinear solver fails
  integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations

  !======= Internals ============

  retval = FARKodeGetNumSteps(arkode_mem, nsteps)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumSteps, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumStepAttempts, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumRhsEvals(arkode_mem, 0, nfe)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumRhsEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumRhsEvals(arkode_mem, 1, nfi)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumRhsEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetActualInitStep(arkode_mem, hinused)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetActualInitStep, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetLastStep(arkode_mem, hlast)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetLastStep, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetCurrentStep(arkode_mem, hcur)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetCurrentStep, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetCurrentTime(arkode_mem, tcur)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetCurrentTime, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumLinSolvSetups, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumErrTestFails(arkode_mem, netfails)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumErrTestFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumNonlinSolvIters(arkode_mem, nniters)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FARKodeGetNumJacEvals(arkode_mem, njacevals)
  if (retval /= 0) then
    print *, 'Error in FARKodeGetNumJacEvals, retval = ', retval, '; halting'
    stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)', 'Total internal steps taken    =', nsteps
  print '(4x,A,i9)', 'Total internal steps attempts =', nst_a
  print '(4x,A,i9)', 'Total rhs exp function calls  =', nfe
  print '(4x,A,i9)', 'Total rhs imp function calls  =', nfi
  print '(4x,A,i9)', 'Total Jacobian function calls =', njacevals
  print '(4x,A,i9)', 'Num lin solver setup calls    =', nlinsetups
  print '(4x,A,i9)', 'Num error test failures       =', netfails
  print '(4x,A,es12.5)', 'First internal step size      =', hinused
  print '(4x,A,es12.5)', 'Last internal step size       =', hlast
  print '(4x,A,es12.5)', 'Next internal step size       =', hcur
  print '(4x,A,es12.5)', 'Current internal time         =', tcur
  print '(4x,A,i9)', 'Num nonlinear solver iters    =', nniters
  print '(4x,A,i9)', 'Num nonlinear solver fails    =', nncfails
  print *, ' '

  return

end subroutine PrintFinalStats
! ----------------------------------------------------------------
