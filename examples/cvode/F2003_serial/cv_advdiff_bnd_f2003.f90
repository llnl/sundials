! ------------------------------------------------------------------
! Programmer(s): Daniel M. Margolis @ UMBC
!                based off the previous Fortran-77 example program,
!                cvode/fcmix_serial/fcvAdvDiff_bnd.f
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
! The following is a simple example problem with a banded
! Jacobian. The problem is the semi-discrete form of the
! advection-diffusion equation in 2D:
! du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
! on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
! interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
! are posed, and the initial condition is the following:
!
!     u(x, y, t = 0) = x*(2 - x)*y*(1 - y)*exp(5*x*y) .
!
! The PDE is discretized on a uniform MX+2 by MY+2 grid with
! central differencing, and with boundary values eliminated,
! leaving an ODE system of size NEQ = MX*MY.
! This program solves this problem with CVODE, using the Fortran/C
! interface routine package. This solution uses the BDF method,
! a user-supplied banded Jacobian routine, and scalar relative and
! absolute tolerances. It prints results at t = .1, .2, ..., 1.0.
! At the end of the run, various counters of interest are printed.
! ------------------------------------------------------------------

module advdiff_mod

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

  ! setup and number of equations
  integer(kind=myindextype), parameter :: mx = 10, my = 5
  integer(kind=myindextype), parameter :: mxmy = mx*my
  integer(kind=myindextype), parameter :: neq = mxmy

  ! ODE constant parameters
  real(c_double), parameter :: xmax = 2.0d0, ymax = 1.0d0
  real(c_double), parameter :: dtout = 0.1d0
  real(c_double), parameter :: dx = xmax/(mx + 1)
  real(c_double), parameter :: dy = ymax/(my + 1)
  real(c_double), parameter :: hdcoef = 1.0d0/(dx*dx)
  real(c_double), parameter :: hacoef = 0.5d0/(2.0d0*dx)
  real(c_double), parameter :: vdcoef = 1.0d0/(dy*dy)

  ! Solving assistance fixed parameters
  real(c_double), parameter :: rtol = 0.0d0
  real(c_double), parameter :: atol = 1.0d-5

  ! ODE non-constant parameters
  integer(kind=myindextype) :: i, j   ! index variables
  integer(kind=myindextype) :: mu, ml ! band preconditioner constants
  real(c_double) :: x, y      ! initialization index variables
  real(c_double) :: unorm     ! solution output variable

contains

  ! ----------------------------------------------------------------
  ! RhsFn provides the right hand side implicit function for the ODE.
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function RhsFn(tn, sunvec_u, sunvec_f, user_data) &
    result(ierr) bind(C, name='RhsFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn         ! current time
    type(N_Vector)        :: sunvec_u   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    type(c_ptr), value :: user_data  ! user-defined data

    ! local data
    real(c_double) :: uij, udn, uup, ult, urt, hdiff, hadv, vdiff

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(mx, my) :: uvec(:, :)
    real(c_double), pointer, dimension(mx, my) :: fvec(:, :)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    uvec(1:mx, 1:my) => FN_VGetArrayPointer(sunvec_u)
    fvec(1:mx, 1:my) => FN_VGetArrayPointer(sunvec_f)

    ! Loop over all grid points
    do i = 1, mx
      do j = 1, my

        ! Extract u at x_i, y_j and four neighboring points.
        uij = uvec(i, j)
        udn = 0.0d0
        if (j /= 1) udn = uvec(i, j - 1)
        uup = 0.0d0
        if (j /= my) uup = uvec(i, j + 1)
        ult = 0.0d0
        if (i /= 1) ult = uvec(i - 1, j)
        urt = 0.0d0
        if (i /= mx) urt = uvec(i + 1, j)

        ! Set diffusion and advection terms and load into fvec.
        hdiff = hdcoef*(ult - 2.0d0*uij + urt)
        hadv = hacoef*(urt - ult)
        vdiff = vdcoef*(uup - 2.0d0*uij + udn)
        fvec(i, j) = hdiff + hadv + vdiff

      end do
    end do

    ! return success
    ierr = 0
    return

  end function RhsFn
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! JacFn provides the user-supplied banded Jacobian
  ! function for the ODE.
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function JacFn(t, sunvec_u, sunvec_f, sunmat_J, &
                                user_data, sunvec_t1, sunvec_t2, sunvec_t3) result(ierr) &
    bind(C, name='JacFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsunmatrix_band_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t
    type(N_Vector)        :: sunvec_u
    type(N_Vector)        :: sunvec_f
    type(SUNMatrix)       :: sunmat_J
    type(c_ptr), value :: user_data
    type(N_Vector)        :: sunvec_t1
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! local data
    integer(kind=myindextype) :: mband, k, ioff, mu1, mu2, smu, mdim
    integer(kind=myindextype) :: start
    real(c_double), pointer, dimension(mdim, neq) :: Jmat(:, :)

    smu = FSUNBandMatrix_StoredUpperBandwidth(sunmat_J)
    mdim = smu + 1 + ml
    Jmat(1:mdim, 1:neq) => FSUNBandMatrix_Data(sunmat_J)

    mu1 = smu + 1
    mu2 = smu + 2
    mband = smu + 1 + ml
    start = smu - mu + 1

    ! Loop over all grid points
    do i = 1, mx
      ioff = (i - 1)*my
      do j = 1, my
        k = j + ioff

        ! Set Jacobian elements in column k of J.
        Jmat(mu1, k) = -2.0d0*(vdcoef + hdcoef)
        if (i /= 1) Jmat(start, k) = hdcoef + hacoef
        if (i /= mx) Jmat(mband, k) = hdcoef - hacoef
        if (j /= 1) Jmat(smu, k) = vdcoef
        if (j /= my) Jmat(mu2, k) = vdcoef

      end do
    end do

    ! return success
    ierr = 0
    return

  end function JacFn
  ! ----------------------------------------------------------------

end module advdiff_mod
! ------------------------------------------------------------------

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fcvode_mod                 ! Fortran interface to the CVode module
  use fsunmatrix_band_mod        ! Fortran interface to banded SUNMatrix
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunlinsol_band_mod        ! Fortran interface to banded SUNLinearSolver
  use advdiff_mod                ! Advection Diffusion functions

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)     :: ctx        ! SUNDIALS context for the simulation
  real(c_double)  :: tstart     ! initial time
  real(c_double)  :: tout       ! output time
  real(c_double)  :: tcur(1)    ! current time
  integer(c_int)  :: ierr       ! error flag from C functions
  integer(c_long) :: outstep    ! output step

  type(N_Vector), pointer :: sunvec_u      ! sundials vector
  type(SUNLinearSolver), pointer :: sunls         ! sundials linear solver
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix (empty)
  type(c_ptr)                    :: cvode_mem     ! CVODE memory
  real(c_double), pointer, dimension(mx, my) :: u(:, :) ! underlying vector

  ! output statistic variables
  integer(c_long)  :: lnst(1)

  !======= Internals ============

  ! create the SUNDIALS context
  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)

  ! initialize ODE
  tstart = 0.0d0
  tcur = tstart
  mu = my
  ml = my

  ! create SUNDIALS N_Vector
  sunvec_u => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_u)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  u(1:mx, 1:my) => FN_VGetArrayPointer(sunvec_u)

  ! initialize and fill initial condition vector
  do i = 1, mx
    x = i*dx
    do j = 1, my
      y = j*dy
      u(i, j) = x*(xmax - x)*y*(ymax - y)*exp(5.0d0*x*y)
    end do
  end do

  ! create and initialize CVode memory
  cvode_mem = FCVodeCreate(CV_BDF, ctx)
  if (.not. c_associated(cvode_mem)) print *, 'ERROR: cvode_mem = NULL'

  ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_u)
  if (ierr /= 0) then
    print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
    stop 1
  end if

  ! Tell CVODE to use a Band linear solver.
  sunmat_A => FSUNBandMatrix(neq, mu, ml, ctx)
  if (.not. associated(sunmat_A)) then
    print *, 'ERROR: sunmat = NULL'
    stop 1
  end if
  sunls => FSUNLinSol_Band(sunvec_u, sunmat_A, ctx)
  if (.not. associated(sunls)) then
    print *, 'ERROR: sunls = NULL'
    stop 1
  end if

  ! Attach the linear solver (with NULL SUNMatrix object)
  ierr = FCVodeSetLinearSolver(cvode_mem, sunls, sunmat_A)
  if (ierr /= 0) then
    print *, 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
  if (ierr /= 0) then
    print *, 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
  if (ierr /= 0) then
    print *, 'Error in FCVodeSetJacFn, ierr = ', ierr, '; halting'
    stop 1
  end if

  ! Start time stepping
  print *, '   '
  print *, 'Band example problem:'
  print '(a,i2)', ' Advection-diffusion, NEQ = ', neq
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '    t       max.norm(u) | lnst'
  print *, ' ------------------------------'

  unorm = maxval(abs(u))
  print '(2x,f6.2,2x,es14.6,2x,i5)', tcur, unorm, lnst

  tout = dtout
  do outstep = 1, 10

    ! call CVode
    ierr = FCVode(cvode_mem, tout, sunvec_u, tcur, CV_NORMAL)
    if (ierr /= 0) then
      print *, 'Error in FCVodeEvolve, ierr = ', ierr, '; halting'
      stop 1
    end if

    ierr = FCVodeGetNumSteps(cvode_mem, lnst)
    if (ierr /= 0) then
      print *, 'Error in FCVodeGetNumSteps, ierr = ', ierr, '; halting'
      stop 1
    end if

    ! print current solution and output statistics
    unorm = maxval(abs(u))
    print '(2x,f6.2,2x,es14.6,2x,i5)', tcur, unorm, lnst

    ! update tout
    tout = tout + dtout

  end do
  print *, ' ------------------------------'

  ! diagnostics output
  call CVodeStats(cvode_mem)

  ! clean up
  call FCVodeFree(cvode_mem)
  call FN_VDestroy(sunvec_u)
  call FSUNMatDestroy(sunmat_A)
  ierr = FSUNLinSolFree(sunls)
  ierr = FSUNContext_Free(ctx)

end program main
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! CVodeStats
!
! Print CVODE statstics to stdandard out
! ----------------------------------------------------------------
subroutine CVodeStats(cvode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use fcvode_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

  integer(c_int)  :: ierr          ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nfe(1)        ! num function evals
  integer(c_long) :: netfails(1)   ! num error test fails
  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nliters(1)    ! linear solver iterations
  integer(c_long) :: ncf(1)        ! num convergence failures nonlinear
  integer(c_long) :: ncfl(1)       ! num convergence failures linear

  !======= Internals ============

  ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumSteps, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeGetNumRhsEvals(cvode_mem, nfe)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumRhsEvals, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeGetNumErrTestFails(cvode_mem, netfails)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumErrTestFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeGetNumNonlinSolvIters(cvode_mem, nniters)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumNonlinSolvIters, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeGetNumLinIters(cvode_mem, nliters)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumLinIters, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeGetNumLinConvFails(cvode_mem, ncfl)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumLinConvFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FCVodeGetNumNonlinSolvConvFails(cvode_mem, ncf)
  if (ierr /= 0) then
    print *, 'Error in FCVodeGetNumNonlinSolvConvFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)', 'Total internal steps taken      =', nsteps
  print '(4x,A,i9)', 'Total rhs function call         =', nfe
  print '(4x,A,i9)', 'Num error test failures         =', netfails
  print '(4x,A,i9)', 'Num nonlinear solver iters      =', nniters
  print '(4x,A,i9)', 'Num linear solver iters         =', nliters
  print '(4x,A,i9)', 'Num nonlinear solver fails      =', ncf
  print '(4x,A,i9)', 'Num linear solver fails         =', ncfl
  print *, ' '

  return

end subroutine CVodeStats
! ----------------------------------------------------------------
