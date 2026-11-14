! ------------------------------------------------------------------
! Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
!                modified by Jean M. Sexton @ LBL
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
! The following is a simple example problem with an analytical
! solution.
!
!   dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
!
! for t in the interval [0.0, 10.0], with initial condition: y=0.
!
! The stiffness of the problem is directly proportional to the
! value of lambda. The value of lambda should be negative to
! result in a well-posed ODE; for values with magnitude larger
! than 100 the problem becomes quite stiff.
! ------------------------------------------------------------------

module analytic_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  !======= Declarations =========
  implicit none

  ! Since SUNDIALS can be compiled with 32-bit or 64-bit sunindextype
  ! we set the integer kind used for indices in this example based
  ! on the the index size SUNDIALS was compiled with so that it works
  ! in both configurations. This is not a requirement for user codes.
#if defined(SUNDIALS_INT32_T)
  integer, parameter :: myindextype = selected_int_kind(8)
#elif defined(SUNDIALS_INT64_T)
  integer, parameter :: myindextype = selected_int_kind(16)
#endif

  ! number of equations
  integer(kind=myindextype), parameter :: neq = 1

  ! ODE parameters
  real(c_double), parameter :: lambda = -100.0d0

contains

  ! ----------------------------------------------------------------
  ! RhsFn provides the right hand side function for the
  ! ODE: dy/dt = f(t,y)
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
    result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn         ! current time
    type(N_Vector)        :: sunvec_y   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    type(c_ptr), value :: user_data  ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yvec(:)
    real(c_double), pointer, dimension(neq) :: fvec(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)

    ! fill RHS vector
    fvec(1) = lambda*yvec(1) + 1.0/(1.0 + tn*tn) - lambda*atan(tn); 
    ! return success
    ierr = 0
    return

  end function RhsFn
  ! ----------------------------------------------------------------

end module analytic_mod
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                ! Fortran interface to the ARKODE
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep time-stepper module
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use fsunadaptcontroller_soderlind_mod ! Fortran interface to Soderlind controller
  use analytic_mod               ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)    :: ctx                      ! SUNDIALS context for the simulation
  real(c_double) :: tstart                   ! initial time
  real(c_double) :: tend                     ! final time
  real(c_double) :: rtol, atol               ! relative and absolute tolerance
  real(c_double) :: dtout                    ! output time interval
  real(c_double) :: tout                     ! output time
  real(c_double) :: tcur(1)                  ! current time
  integer(c_int) :: ierr                     ! error flag from C functions
  integer(c_int) :: nout                     ! number of outputs
  integer(c_int) :: outstep                  ! output loop counter

  type(N_Vector), pointer :: sunvec_y    ! sundials vector
  type(SUNMatrix), pointer :: sunmat_A    ! sundials matrix
  type(SUNLinearSolver), pointer :: sunls       ! sundials linear solver
  type(SUNAdaptController), pointer :: sunCtrl     ! time step controller
  type(c_ptr)                       :: arkode_mem  ! ARKODE memory
  real(c_double), pointer, dimension(neq) :: yvec(:)  ! underlying vector

  !======= Internals ============

  ! create the SUNDIALS context
  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)

  ! initialize ODE
  tstart = 0.0d0
  tend = 10.0d0
  tcur = tstart
  tout = tstart
  dtout = 1.0d0
  nout = ceiling(tend/dtout)

  ! create SUNDIALS N_Vector
  sunvec_y => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_y)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  yvec => FN_VGetArrayPointer(sunvec_y)

  ! initialize solution vector
  call FN_VConst(0.0d0, sunvec_y)

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(RhsFn), tstart, sunvec_y, ctx)
  if (.not. c_associated(arkode_mem)) print *, 'ERROR: arkode_mem = NULL'

  ! Tell ARKODE to use a dense linear solver.
  sunmat_A => FSUNDenseMatrix(neq, neq, ctx)
  if (.not. associated(sunmat_A)) then
    print *, 'ERROR: sunmat = NULL'
    stop 1
  end if

  sunls => FSUNLinSol_Dense(sunvec_y, sunmat_A, ctx)
  if (.not. associated(sunls)) then
    print *, 'ERROR: sunls = NULL'
    stop 1
  end if

  ierr = FARKodeSetLinearSolver(arkode_mem, sunls, sunmat_A)
  if (ierr /= 0) then
    write (*, *) 'Error in FARKodeSetLinearSolver'
    stop 1
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-6
  atol = 1.0d-10

  ierr = FARKodeSStolerances(arkode_mem, rtol, atol)
  if (ierr /= 0) then
    write (*, *) 'Error in FARKodeSStolerances, ierr = ', ierr, '; halting'
    stop 1
  end if

  sunCtrl => FSUNAdaptController_ImpGus(ctx)
  if (.not. associated(sunCtrl)) then
    print *, 'ERROR: sunCtrl = NULL'
    stop 1
  end if
  ierr = FARKodeSetAdaptController(arkode_mem, sunCtrl)
  if (ierr /= 0) then
    write (*, *) 'Error in FARKodeSetAdaptController, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeSetNonlinConvCoef(arkode_mem, 0.01d0)
  if (ierr /= 0) then
    write (*, *) 'Error in FARKodeSetNonlinConvCoef, ierr = ', ierr, '; halting'
    stop 1
  end if

  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '       t           y        '
  print *, '----------------------------'
  print '(2x,2(es12.5,1x))', tcur, yvec(1)
  do outstep = 1, nout

    ! call ARKodeEvolve
    tout = min(tout + dtout, tend)
    ierr = FARKodeEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
    if (ierr /= 0) then
      write (*, *) 'Error in FARKODE, ierr = ', ierr, '; halting'
      stop 1
    end if

    ! output current solution
    print '(2x,2(es12.5,1x))', tcur, yvec(1)

  end do

  ! diagnostics output
  call ARKStepStats(arkode_mem)

  ! clean up
  call FARKodeFree(arkode_mem)
  call FN_VDestroy(sunvec_y)
  call FSUNMatDestroy(sunmat_A)
  ierr = FSUNLinSolFree(sunls)
  ierr = FSUNAdaptController_Destroy(sunCtrl)
  ierr = FSUNContext_Free(ctx)

end program main
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! ARKStepStats
!
! Print ARKODE statstics to stdandard out
! ----------------------------------------------------------------
subroutine ARKStepStats(arkode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: arkode_mem ! solver memory structure

  integer(c_int)  :: ierr          ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nst_a(1)      ! num steps attempted
  integer(c_long) :: nfevals(1)    ! num function evals
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

  ierr = FARKodeGetNumSteps(arkode_mem, nsteps)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumSteps, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumStepAttempts, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumRhsEvals(arkode_mem, -1, nfevals)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumRhsEvals, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetActualInitStep(arkode_mem, hinused)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetActualInitStep, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetLastStep(arkode_mem, hlast)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetLastStep, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetCurrentStep(arkode_mem, hcur)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetCurrentStep, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetCurrentTime(arkode_mem, tcur)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetCurrentTime, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumLinSolvSetups, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumErrTestFails(arkode_mem, netfails)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumErrTestFails, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumNonlinSolvIters(arkode_mem, nniters)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvIters, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvConvFails, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumJacEvals(arkode_mem, njacevals)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumJacEvals, retval = ', ierr, '; halting'
    stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)', 'Total internal steps taken    =', nsteps
  print '(4x,A,i9)', 'Total internal steps attempts =', nst_a
  print '(4x,A,i9)', 'Total rhs function calls      =', nfevals
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

end subroutine ARKStepStats
! ----------------------------------------------------------------
