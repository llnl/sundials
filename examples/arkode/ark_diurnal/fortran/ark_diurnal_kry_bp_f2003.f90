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
! FARKODE Example Problem: 2D kinetics-transport,
! precond. Krylov solver.
!
! An ODE system is generated from the following 2-species diurnal
! kinetics advection-diffusion PDE system in 2 space dimensions:
!
! dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
!                       + Ri(c1,c2,t)      for i = 1,2,   where
! R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
! R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
! Kv(y) = Kv0*exp(y/5) ,
! Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
! vary diurnally.
!
! The problem is posed on the square
! 0 <= x <= 20,    30 <= y <= 50   (all in km),
! with homogeneous Neumann boundary conditions, and for time t in
! 0 <= t <= 86400 sec (1 day).
!
! The PDE system is treated by central differences on a uniform
! 10 x 10 mesh, with simple polynomial initial profiles.
! The problem is solved with ARKODE, with the DIRK/GMRES method and
! using the FARKBP banded preconditioner.
!
! Note that this problem should only work with SUNDIALS configured
! to use 'realtype' as 'double' and 'sunindextype' as '64bit'
!
! The second and third dimensions of U here must match the values of
! MX and MY, for consistency with the output statements below.
! ------------------------------------------------------------------

module DiurnalKryBP_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  !======= Declarations =========
  implicit none

  ! setup and number of equations
  integer(c_int), parameter :: mx = 10, my = 10
  integer(c_int64_t), parameter :: mm = mx*my
  integer(c_int64_t), parameter :: neq = 2*mm

  ! ODE constant parameters
  real(c_double), parameter :: Kh = 4.0d-6
  real(c_double), parameter :: Vel = 0.001d0
  real(c_double), parameter :: Kv0 = 1.0d-8
  real(c_double), parameter :: q1 = 1.63d-16
  real(c_double), parameter :: q2 = 4.66d-16
  real(c_double), parameter :: c3 = 3.7d16
  real(c_double), parameter :: pi = 3.1415926535898d0
  real(c_double), parameter :: halft = 4.32d4
  real(c_double), parameter :: om = pi/halft
  real(c_double), parameter :: dx = 20.0d0/(mx - 1.0d0)
  real(c_double), parameter :: dy = 20.0d0/(my - 1.0d0)
  real(c_double), parameter :: hdco = Kh/(dx**2)
  real(c_double), parameter :: haco = Vel/(2.0d0*dx)
  real(c_double), parameter :: vdco = (1.0d0/(dy**2))*Kv0
  real(c_double), parameter :: a3 = 22.62d0
  real(c_double), parameter :: a4 = 7.601d0

  ! Solving assistance fixed parameters
  real(c_double), parameter  :: twohr = 7200.0D0
  real(c_double), parameter  :: rtol = 1.0d-5
  real(c_double), parameter  :: floor = 100.0d0
  real(c_double), parameter  :: delta = 0.0d0
  real(c_double), parameter  :: atol = rtol*floor
  integer(c_int), parameter  :: Jpretype = 1
  integer(c_int), parameter  :: iGStype = 1
  integer(c_int), parameter  :: maxL = 0
  integer(c_long), parameter :: mxsteps = 10000

  ! ODE non-constant parameters
  real(c_double) :: q3
  real(c_double) :: q4
  real(c_double) :: c1
  real(c_double) :: c2
  integer(c_int) :: jx, jy

contains

  ! ----------------------------------------------------------------
  ! ImpRhsFn provides the right hand side implicit function for the
  ! ODE: dy1/dt = f1(t,y1,y2,y3)
  !      dy2/dt = f2(t,y1,y2,y3)
  !      dy3/dt = f3(t,y1,y2,y3)
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function ImpRhsFn(tn, sunvec_u, sunvec_f, user_data) &
    result(ierr) bind(C, name='ImpRhsFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_u  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! local data
    integer(c_int)  :: jleft, jright, jdn, jup
    real(c_double)  :: c1dn, c2dn, c1up, c2up, c1lt, c2lt
    real(c_double)  :: c1rt, c2rt, cydn, cyup, hord1, hord2, horad1
    real(c_double)  :: horad2, qq1, qq2, qq3, qq4, rkin1, rkin2, s
    real(c_double)  :: vertd1, vertd2, ydn, yup

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(2, mx, my) :: uvecI(:, :, :)
    real(c_double), pointer, dimension(2, mx, my) :: fvecI(:, :, :)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    uvecI(1:2, 1:mx, 1:my) => FN_VGetArrayPointer(sunvec_u)
    fvecI(1:2, 1:mx, 1:my) => FN_VGetArrayPointer(sunvec_f)

    ! Set diurnal rate coefficients.
    s = sin(om*tn)
    if (s > 0.0d0) then
      q3 = exp(-a3/s)
      q4 = exp(-a4/s)
    else
      q3 = 0.0d0
      q4 = 0.0d0
    end if

    ! Loop over all grid points.
    do jy = 1, my
      ydn = 30.0d0 + (jy - 1.5d0)*dy
      yup = ydn + dy
      cydn = vdco*exp(0.2d0*ydn)
      cyup = vdco*exp(0.2d0*yup)
      jdn = jy - 1
      if (jy == 1) jdn = my
      jup = jy + 1
      if (jy == my) jup = 1
      do jx = 1, mx
        c1 = uvecI(1, jx, jy)
        c2 = uvecI(2, jx, jy)
        ! Set kinetic rate terms.
        qq1 = q1*c1*c3
        qq2 = q2*c1*c2
        qq3 = q3*c3
        qq4 = q4*c2
        rkin1 = -qq1 - qq2 + 2.0d0*qq3 + qq4
        rkin2 = qq1 - qq2 - qq4
        ! Set vertical diffusion terms.
        c1dn = uvecI(1, jx, jdn)
        c2dn = uvecI(2, jx, jdn)
        c1up = uvecI(1, jx, jup)
        c2up = uvecI(2, jx, jup)
        vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn)
        vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn)
        ! Set horizontal diffusion and advection terms.
        jleft = jx - 1
        if (jx == 1) jleft = mx
        jright = jx + 1
        if (jx == mx) jright = 1
        c1lt = uvecI(1, jleft, jy)
        c2lt = uvecI(2, jleft, jy)
        c1rt = uvecI(1, jright, jy)
        c2rt = uvecI(2, jright, jy)
        hord1 = hdco*(c1rt - 2.0d0*c1 + c1lt)
        hord2 = hdco*(c2rt - 2.0d0*c2 + c2lt)
        horad1 = haco*(c1rt - c1lt)
        horad2 = haco*(c2rt - c2lt)
        ! load all terms into fvecI.
        fvecI(1, jx, jy) = vertd1 + hord1 + horad1 + rkin1
        fvecI(2, jx, jy) = vertd2 + hord2 + horad2 + rkin2
      end do
    end do

    ! return success
    ierr = 0
    return

  end function ImpRhsFn
  ! ----------------------------------------------------------------

end module DiurnalKryBP_mod
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                ! Fortran interface to the ARKODE module
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep module
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunlinsol_spgmr_mod       ! Fortran interface to spgmr SUNLinearSolver
  use DiurnalKryBP_mod           ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)     :: ctx        ! SUNDIALS context for the simulation
  real(c_double)  :: tstart     ! initial time
  real(c_double)  :: tout       ! output time
  real(c_double)  :: tcur(1)    ! current time
  real(c_double)  :: cx, cy     ! initialization variables
  integer(c_int)  :: ierr       ! error flag from C functions
  integer(c_long) :: outstep    ! output step
  integer(c_int64_t) :: mu, ml     ! band preconditioner constants
  real(c_double)  :: x, y       ! initialization index variables

  type(N_Vector), pointer :: sunvec_u      ! sundials vector
  type(N_Vector), pointer :: sunvec_f      ! sundials vector
  type(SUNLinearSolver), pointer :: sunls         ! sundials linear solver
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix (empty)
  type(c_ptr)                    :: arkode_mem    ! ARKODE memory
  real(c_double), pointer, dimension(2, mx, my) :: uvec(:, :, :) ! underlying vector

  ! output statistic variables
  integer(c_long) :: lnst(1), lnst_att(1)
  real(c_double)  :: lh(1)

  !======= Internals ============

  ! create the SUNDIALS context
  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)

  ! initialize ODE
  tstart = 0.0d0
  tcur = tstart

  ! create SUNDIALS N_Vectors
  sunvec_u => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_u)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  uvec(1:2, 1:mx, 1:my) => FN_VGetArrayPointer(sunvec_u)

  sunvec_f => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_f)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  ! initialize and fill initial condition vector
  do jy = 1, my
    y = 30.0d0 + (jy - 1.0d0)*dy
    cy = (0.1d0*(y - 40.0d0))**2
    cy = 1.0d0 - cy + 0.5d0*cy**2
    do jx = 1, mx
      x = (jx - 1.0d0)*dx
      cx = (0.1d0*(x - 10.0d0))**2
      cx = 1.0d0 - cx + 0.5d0*cx**2
      uvec(1, jx, jy) = 1.0d6*cx*cy
      uvec(2, jx, jy) = 1.0d12*cx*cy
    end do
  end do

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(ImpRhsFn), tstart, sunvec_u, ctx)
  if (.not. c_associated(arkode_mem)) print *, 'ERROR: arkode_mem = NULL'

  ! Tell ARKODE to use a SPGMR linear solver.
  sunls => FSUNLinSol_SPGMR(sunvec_u, Jpretype, maxL, ctx)
  if (.not. associated(sunls)) then
    print *, 'ERROR: sunls = NULL'
    stop 1
  end if

  ! Attach the linear solver (with NULL SUNMatrix object)
  sunmat_A => null()
  ierr = FARKodeSetLinearSolver(arkode_mem, sunls, sunmat_A)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetLinearSolver, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FSUNLinSol_SPGMRSetGSType(sunls, iGStype)
  if (ierr /= 0) then
    print *, 'Error in FSUNLinSol_SPGMRSetGSType, ierr = ', ierr, '; halting'
    stop 1
  end if

  mu = 2
  ml = 2
  ierr = FARKBandPrecInit(arkode_mem, neq, mu, ml)
  if (ierr /= 0) then
    print *, 'Error in FARKBandPrecInit, ierr = ', ierr, '; halting'
    stop 1
  end if

  ! Set additional method parameters
  ierr = FARKodeSStolerances(arkode_mem, rtol, atol)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSStolerances, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeSetMaxNumSteps(arkode_mem, mxsteps)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetMaxNumSteps'
    stop 1
  end if

  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '      t          c1  (bottom left      middle       top right)  | lnst  lnst_att  lh'
  print *, '      t          c2  (bottom left      middle       top right)  | lnst  lnst_att  lh'
  print *, ' ----------------------------------------------------------------------------------------'
  tout = twohr
  do outstep = 1, 12

    ! call ARKodeEvolve
    ierr = FARKodeEvolve(arkode_mem, tout, sunvec_u, tcur, ARK_NORMAL)
    if (ierr /= 0) then
      print *, 'Error in FARKodeEvolve, ierr = ', ierr, '; halting'
      stop 1
    end if

    ! get some solver statistics
    ierr = FARKodeGetNumSteps(arkode_mem, lnst)
    if (ierr /= 0) then
      print *, 'Error in FARKodeGetNumSteps, ierr = ', ierr, '; halting'
      stop 1
    end if

    ierr = FARKodeGetNumStepAttempts(arkode_mem, lnst_att)
    if (ierr /= 0) then
      print *, 'Error in FARKodeGetNumStepAttempts, ierr = ', ierr, '; halting'
      stop 1
    end if

    ierr = FARKodeGetCurrentStep(arkode_mem, lh)
    if (ierr /= 0) then
      print *, 'Error in FARKodeGetCurrentStep, ierr = ', ierr, '; halting'
      stop 1
    end if

    ! print current solution and output statistics
    print '(2x,4(es14.6,2x),i5,i5,es14.6)', tcur, uvec(1, 1, 1), uvec(1, 5, 5), uvec(1, 10, 10), lnst, lnst_att, lh
    print '(18x,3(es14.6,2x))', uvec(2, 1, 1), uvec(2, 5, 5), uvec(2, 10, 10)

    ! update tout
    tout = tout + twohr

  end do
  print *, ' ----------------------------------------------------------------------------------------'

  ! diagnostics output
  call ARKStepStats(arkode_mem)

  ! clean up
  call FARKodeFree(arkode_mem)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_f)
  ierr = FSUNLinSolFree(sunls)
  ierr = FSUNContext_Free(ctx)

end program main
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! ARKStepStats
!
! Print ARKODE statstics to standard out
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
  integer(c_long) :: nfe(1)        ! num explicit function evals
  integer(c_long) :: nfi(1)        ! num implicit function evals
  integer(c_long) :: netfails(1)   ! num error test fails
  integer(c_long) :: npe(1)        ! num preconditioner evals
  integer(c_long) :: nps(1)        ! num preconditioner solves
  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nliters(1)    ! linear solver iterations
  integer(c_long) :: ncf(1)        ! num convergence failures nonlinear
  integer(c_long) :: ncfl(1)       ! num convergence failures linear
  integer(c_long) :: lenrw(1)      ! main solver real/int workspace size
  integer(c_long) :: leniw(1)
  integer(c_long) :: lenrwls(1)    ! linear solver real/int workspace size
  integer(c_long) :: leniwls(1)
  integer(c_long) :: nfebp(1)      ! num f evaluations
  real(c_double)  :: avdim(1)      ! avg Krylov subspace dim (NLI/NNI)
  integer(c_long) :: lenrwbp(1)    ! band preconditioner real/int workspace size
  integer(c_long) :: leniwbp(1)

  !======= Internals ============

  ierr = FARKodeGetNumSteps(arkode_mem, nsteps)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumSteps, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumStepAttempts, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumRhsEvals(arkode_mem, 0, nfe)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumRhsEvals, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumRhsEvals(arkode_mem, 1, nfi)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumRhsEvals, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumErrTestFails(arkode_mem, netfails)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumErrTestFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumPrecEvals(arkode_mem, npe)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumPrecEvals, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumPrecSolves(arkode_mem, nps)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumPrecSolves, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumNonlinSolvIters(arkode_mem, nniters)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvIters, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumLinIters(arkode_mem, nliters)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumLinIters, ierr = ', ierr, '; halting'
    stop 1
  end if

  avdim = dble(nliters)/dble(nniters)

  ierr = FARKodeGetNumLinConvFails(arkode_mem, ncfl)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumLinConvFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumNonlinSolvConvFails(arkode_mem, ncf)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvConvFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetWorkSpace(arkode_mem, lenrw, leniw)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetWorkSpace, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetLinWorkSpace(arkode_mem, lenrwls, leniwls)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetLinWorkSpace, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKBandPrecGetWorkSpace(arkode_mem, lenrwbp, leniwbp)
  if (ierr /= 0) then
    print *, 'Error in FARKBandPrecGetWorkSpace, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKBandPrecGetNumRhsEvals(arkode_mem, nfebp)
  if (ierr /= 0) then
    print *, 'Error in FARKBandPrecGetNumRhsEvals, ierr = ', ierr, '; halting'
    stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)', 'Total internal steps taken      =', nsteps
  print '(4x,A,i9)', 'Total internal steps attempts   =', nst_a
  print '(4x,A,i9)', 'Total rhs exp function call     =', nfe
  print '(4x,A,i9)', 'Total rhs imp function call     =', nfi
  print '(4x,A,i9)', 'Total num preconditioner evals  =', npe
  print '(4x,A,i9)', 'Total num preconditioner solves =', nps
  print '(4x,A,i9)', 'Num error test failures         =', netfails
  print '(4x,A,i9)', 'Num nonlinear solver iters      =', nniters
  print '(4x,A,i9)', 'Num linear solver iters         =', nliters
  print '(4x,A,es14.6)', 'Avg Krylov subspace dim         =', avdim
  print '(4x,A,i9)', 'Num nonlinear solver fails      =', ncf
  print '(4x,A,i9)', 'Num linear solver fails         =', ncfl
  print '(4x,A,2(i9,3x))', 'main solver real/int workspace sizes   =', lenrw, leniw
  print '(4x,A,2(i9,3x))', 'linear solver real/int workspace sizes =', lenrwls, leniwls
  print '(4x,A,2(i9,3x))', 'ARKBandPre real/int workspace sizes    =', lenrwbp, leniwbp
  print '(4x,A,i9)', 'ARKBandPre number of f evaluations     =', nfebp
  print *, ' '

  return

end subroutine ARKStepStats
! ----------------------------------------------------------------
