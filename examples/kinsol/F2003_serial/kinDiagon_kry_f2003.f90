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
! Simple diagonal test with Fortran interface, using user-supplied
! preconditioner setup and solve routines (supplied in Fortran).
!
! This example does a basic test of the solver by solving the
! system:
!    f(u) = 0  for
!    f(u) = u(i)^2 - i^2
!
! No scaling is done.
! An approximate diagonal preconditioner is used.
!
! ------------------------------------------------------------------

module kinDiagonKry_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod ! access sundials types

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
  integer(kind=myindextype), parameter :: neq = 128

  ! reciprocal of the Jacobian diagonal
  real(c_double) :: p(neq)

contains

  ! ----------------------------------------------------------------
  ! init: Initializes variables u, scale, constr
  ! ----------------------------------------------------------------
  subroutine init(sunvec_u, sunvec_s, sunvec_c)

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector) :: sunvec_u  ! solution N_Vector
    type(N_Vector) :: sunvec_s  ! scaling N_Vector
    type(N_Vector) :: sunvec_c  ! constraint N_Vector

    ! local variables
    integer(c_int64_t) :: i

    ! pointers to data arrays in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: u(:), scale(:), constr(:)

    u(1:neq) => FN_VGetArrayPointer(sunvec_u)
    scale(1:neq) => FN_VGetArrayPointer(sunvec_s)
    constr(1:neq) => FN_VGetArrayPointer(sunvec_c)

    ! Set initial guess, disable scaling and constraints
    do i = 1, neq
      u(i) = 2.0d0*dble(i)
    end do
    scale = 1.0d0
    constr = 0.0d0

  end subroutine init

  ! ----------------------------------------------------------------
  ! func: The nonlinear residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function func(sunvec_u, sunvec_f, user_data) &
    result(ierr) bind(C)

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_u  ! solution N_Vector
    type(N_Vector)     :: sunvec_f  ! LHS N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! local variables
    integer(c_int64_t) :: i

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: uu(:), ff(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    uu(1:neq) => FN_VGetArrayPointer(sunvec_u)
    ff(1:neq) => FN_VGetArrayPointer(sunvec_f)

    ! loop over domain, computing our system f(u) = 0
    do i = 1, neq
      ff(i) = uu(i)*uu(i) - dble(i*i)
    end do

    ! return success
    ierr = 0
    return

  end function func

  ! ----------------------------------------------------------------
  ! kpsetup: The KINSOL Preconditioner setup function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function kpsetup(sunvec_u, sunvec_s, sunvec_f, &
                                  sunvec_fs, user_data) result(ierr) bind(C)

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_u  ! solution N_Vector
    type(N_Vector)     :: sunvec_s  ! scaling N_Vector
    type(N_Vector)     :: sunvec_f  ! LHS N_Vector
    type(N_Vector)     :: sunvec_fs ! LHS scaling N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! loop counter
    integer(c_int64_t) :: i

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: udata(:)

    !======= Internals ============

    ! get the data array from the SUNDIALS vector
    udata(1:neq) => FN_VGetArrayPointer(sunvec_u)

    ! setup preconditioner
    do i = 1, neq
      p(i) = 0.5d0/(udata(i) + 5.0d0)
    end do

    ! return success
    ierr = 0
    return

  end function kpsetup

  ! ----------------------------------------------------------------
  ! kpsolve: The KINSOL Preconditioner solve function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function kpsolve(sunvec_u, sunvec_s, sunvec_f, &
                                  sunvec_fs, sunvec_v, user_data) result(ierr) bind(C)

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_u  ! solution N_Vector
    type(N_Vector)     :: sunvec_s  ! scaling N_Vector
    type(N_Vector)     :: sunvec_f  ! LHS N_Vector
    type(N_Vector)     :: sunvec_fs ! LHS scaling N_Vector
    type(N_Vector)     :: sunvec_v  ! LHS scaling N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! loop counter
    integer(c_int64_t) :: i

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: v(:)

    !======= Internals ============
    ! get data arrays from SUNDIALS vectors
    v(1:neq) => FN_VGetArrayPointer(sunvec_v)

    ! preconditioner solver
    do i = 1, neq
      v(i) = v(i)*p(i)
    end do

    ! return success
    ierr = 0
    return

  end function kpsolve

end module kinDiagonKry_mod

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use fkinsol_mod                ! Fortran interface to KINSOL
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunlinsol_spgmr_mod       ! Fortran interface to SPGMR SUNLinearSolver
  use kinDiagonKry_mod           ! problem-defining functions

  !======= Declarations =========
  implicit none

  integer(c_int), parameter :: prectype = SUN_PREC_RIGHT ! right precondition
  integer(c_int), parameter :: maxl = 10                 ! Krylov space size
  integer(c_int), parameter :: maxlrst = 2               ! number of restarts
  integer(c_long), parameter :: msbpre = 5               ! iters between setups
  real(c_double), parameter :: fnormtol = 1.0d-5         ! residual tol
  real(c_double), parameter :: scsteptol = 1.0d-4        ! scaled step tol

  ! local variables
  integer(c_int)          :: retval
  type(c_ptr)             :: sunctx               ! sundials context
  type(N_Vector), pointer :: sunvec_u             ! solution vector
  type(N_Vector), pointer :: sunvec_s             ! scaling vector
  type(N_Vector), pointer :: sunvec_c             ! constraint vector
  type(SUNMatrix), pointer :: sunmat_J            ! Jacobian matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! linear solver
  type(c_ptr) :: kmem                             ! KINSOL memory

  !======= Internals ============

  ! -------------------------
  ! Print problem description

  print '(a)', "Example program kinDiagon_kry_f2003:"
  print '(a)', ""
  print '(a)', "This example demonstrates using the KINSOL Fortran interface to"
  print '(a)', "solve a diagonal algebraic system using a Newton-Krylov method"
  print '(a)', "with a diagonal preconditioner in a serial environment."
  print '(a)', ""
  print '(a,i3)', "Problem size: ", neq

  ! -------------------------
  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
  if (retval /= 0) then
    print *, 'Error in FSUNContext_Create'
    stop 1
  end if

  ! -------------------------
  ! Create vectors for solution and scales

  sunvec_u => FN_VNew_Serial(neq, sunctx)
  sunvec_s => FN_VNew_Serial(neq, sunctx)
  sunvec_c => FN_VNew_Serial(neq, sunctx)

  call init(sunvec_u, sunvec_s, sunvec_c)

  ! -------------------------
  ! Initialize and allocate memory for KINSOL

  kmem = FKINCreate(sunctx)
  if (.not. c_associated(kmem)) then
    print *, 'Error kmem = NULL'
    stop 1
  end if

  ! sunvec_u is used as a template

  retval = FKINInit(kmem, c_funloc(func), sunvec_u)
  if (retval /= 0) then
    print *, 'Error in FKINInit, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Set optional inputs

  retval = FKINSetMaxSetupCalls(kmem, msbpre)
  if (retval /= 0) then
    print *, 'Error in FKINSetMaxSetupCalls, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINSetFuncNormTol(kmem, fnormtol)
  if (retval /= 0) then
    print *, 'Error in FKINSetFuncNormTol, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINSetScaledStepTol(kmem, scsteptol)
  if (retval /= 0) then
    print *, 'Error in FKINSetScaledStepTol, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINSetConstraints(kmem, sunvec_c)
  if (retval /= 0) then
    print *, 'Error in FKINSetConstraints, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Create a SPGMR linear solver

  sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_u, prectype, maxl, sunctx)
  if (.not. associated(sunlinsol_LS)) then
    print *, 'Error sunlinsol = NULL'
    stop 1
  end if

  ! -------------------------
  ! Attach linear solver

  sunmat_J => null()

  retval = FKINSetLinearSolver(kmem, sunlinsol_LS, sunmat_J)
  if (retval /= 0) then
    print *, 'Error in FKINSetLinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Set more optional SPGMR inputs

  retval = FSUNLinSol_SPGMRSetMaxRestarts(sunlinsol_LS, maxlrst)
  if (retval /= 0) then
    print *, 'Error in FSUNLinSol_SPGMRSetMaxRestarts, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Set preconditioner functions

  retval = FKINSetPreconditioner(kmem, c_funloc(kpsetup), c_funloc(kpsolve))
  if (retval /= 0) then
    print *, 'Error in FKINSetPreconditioner, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Call KINSol to solve problem
  !
  ! arguments: KINSol memory block
  !            Initial guess on input, solution on output
  !            Globalization strategy choice
  !            Scaling vector for the solution
  !            Scaling vector for the residual

  retval = FKINSol(kmem, sunvec_u, KIN_NONE, sunvec_s, sunvec_s)
  if (retval /= 0) then
    print *, 'Error in FKINSol, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Print solution and solver statistics

  print '(a)', ""
  call PrintOutput(sunvec_u)
  call PrintFinalStats(kmem)

  ! clean up
  call FKINFree(kmem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_s)
  call FN_VDestroy(sunvec_c)
  retval = FSUNContext_Free(sunctx)

end program main

! ----------------------------------------------------------------
! Print solution at selected points
! ----------------------------------------------------------------
subroutine PrintOutput(sunvec_u)

  !======= Inclusions ===========
  use kinDiagonKry_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector) :: sunvec_u

  ! local variables
  real(c_double), pointer, dimension(neq) :: uu(:)
  integer(c_int64_t) :: i

  !======= Internals ============

  uu(1:neq) => FN_VGetArrayPointer(sunvec_u)

  do i = 1, neq, 4
    print '(i4, 4(1x, f10.6))', i, uu(i), uu(i + 1), uu(i + 2), uu(i + 3)
  end do

  return

end subroutine PrintOutput

! ----------------------------------------------------------------
! Print KINSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(kmemo)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fkinsol_mod

  !======= Declarations =========
  implicit none

  type(c_ptr) :: kmemo

  integer(c_int) :: retval
  integer(c_long) :: nni(1), nli(1), nfe(1), npe(1), nps(1), ncfl(1)

  !======= Internals ============

  ! Main solver statistics

  retval = FKINGetNumNonlinSolvIters(kmemo, nni)
  if (retval /= 0) then
    print *, 'Error in FKINGetNumNonlinSolvIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumLinIters(kmemo, nli)
  if (retval /= 0) then
    print *, 'Error in FKINGetNumLinIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumFuncEvals(kmemo, nfe)
  if (retval /= 0) then
    print *, 'Error in FKINGetNumFuncEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumPrecEvals(kmemo, npe)
  if (retval /= 0) then
    print *, 'Error in KINGetNumPrecEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumPrecSolves(kmemo, nps)
  if (retval /= 0) then
    print *, 'Error in KINGetNumPrecSolves, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumLinConvFails(kmemo, ncfl)
  if (retval /= 0) then
    print *, 'Error in KINGetNumLinConvFails, retval = ', retval, '; halting'
    stop 1
  end if

  print '(a)', ''
  print '(a)', 'Final Statistics..'
  print '(a)', ''
  print '(2(A,i6))', 'nni      =', nni, '    nli     =', nli
  print '(2(A,i6))', 'nfe      =', nfe, '    npe     =', npe
  print '(2(A,i6))', 'nps      =', nps, '    nlcf    =', ncfl

  return

end subroutine PrintFinalStats
