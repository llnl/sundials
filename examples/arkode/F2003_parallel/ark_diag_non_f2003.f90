!-----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ UMBC
!                modified by Daniel M. Margolis @ UMBC
!-----------------------------------------------------------------
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
!-----------------------------------------------------------------
! Example problem:
!
! Diagonal ODE example. Nonstiff case: alpha = 10/NEQ.
! Uses ERKStep to explicitly solve nonstiff ODE.
!
!-----------------------------------------------------------------

module DiagnonData
  !---------------------------------------------------------------
  ! Description:
  !    Module containing problem-defining parameters.
  !---------------------------------------------------------------

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  save

  ! SUNDIALS simulation context
  type(c_ptr) :: sunctx

  ! MPI domain decomposition information
  integer, target :: comm  ! communicator object
  integer :: myid          ! MPI process ID
  integer :: nprocs        ! total number of MPI processes

  ! Problem parameters
  integer(c_int64_t), parameter :: nlocal = 2
  integer(c_int64_t)  :: neq
  real(c_double) :: alpha

contains

  !-----------------------------------------------------------------
  ! ODE RHS function f(t,y).
  !-----------------------------------------------------------------
  integer(c_int) function frhs(t, sunvec_y, sunvec_ydot, user_data) &
    result(retval) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t            ! current time
    type(N_Vector)        :: sunvec_y     ! solution N_Vector
    type(N_Vector)        :: sunvec_ydot  ! rhs N_Vector
    type(c_ptr), value :: user_data    ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nlocal) :: y(:)
    real(c_double), pointer, dimension(nlocal) :: ydot(:)

    ! local data
    integer :: i

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    y(1:nlocal) => FN_VGetArrayPointer(sunvec_y)
    ydot(1:nlocal) => FN_VGetArrayPointer(sunvec_ydot)

    ! Initialize ydot to zero
    ydot = 0.d0

    ! Fill ydot with rhs function
    do i = 1, nlocal
      ydot(i) = -alpha*(myid*nlocal + i)*y(i)
    end do

    retval = 0              ! Return with success
    return
  end function frhs
  !-----------------------------------------------------------------

end module DiagnonData
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! Main driver program
!-------------------------------------------------------------------
program driver

  ! inclusions
  use, intrinsic :: iso_c_binding
  use farkode_mod                ! Access ARKODE
  use farkode_erkstep_mod        ! Access ERKStep
  use fnvector_parallel_mod      ! Access parallel N_Vector
  use DiagnonData

  !======= Declarations =========
  implicit none

  ! Declarations
  ! general problem parameters
  integer, parameter :: Nt = 10                 ! total number of output times
  real(c_double), parameter :: T0 = 0.d0        ! initial time
  real(c_double), parameter :: Tf = 1.d0        ! final time
  real(c_double), parameter :: rtol = 1.d-5     ! relative and absolute tolerances
  real(c_double), parameter :: atol = 1.d-10

  ! solution vector and other local variables
  type(N_Vector), pointer :: sunvec_y                ! solution N_Vector
  real(c_double), pointer, dimension(nlocal) :: y(:) ! vector data
  type(c_ptr)     :: arkode_mem                      ! ARKODE memory
  integer(c_int) :: retval
  integer :: ierr
  logical :: outproc
  real(c_double) :: t(1), dTout, tout
  integer(c_long) :: nst(1)     ! number of time steps
  integer(c_long) :: nst_a(1)   ! number of step attempts
  integer(c_long) :: nfe(1)     ! number of explicit RHS evals
  integer(c_long) :: netf(1)    ! number of error test fails
  integer :: i, ioutput
  real(c_double) :: errmax, erri, gerrmax

  ! Initialize MPI variables
  comm = MPI_COMM_WORLD
  myid = 0
  nprocs = 0

  ! initialize MPI
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Init = ", ierr
    stop 1
  end if
  call MPI_Comm_size(comm, nprocs, ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Comm_size = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if
  call MPI_Comm_rank(comm, myid, ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Comm_rank = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Set input arguments neq and alpha
  neq = nprocs*nlocal
  alpha = 10.0d0/neq

  ! Create SUNDIALS simulation context, now that comm has been configured
  retval = FSUNContext_Create(comm, sunctx)
  if (retval /= 0) then
    print *, "Error: FSUNContext_Create returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Initial problem output
  outproc = (myid == 0)
  if (outproc) then
    write (6, *) "  "
    write (6, *) "Diagonal test problem:"; 
    write (6, '(A,i4)') "   neq = ", neq
    write (6, '(A,i4)') "   nlocal = ", nlocal
    write (6, '(A,i4)') "   nprocs = ", nprocs
    write (6, '(A,es9.2)') "   rtol = ", rtol
    write (6, '(A,es9.2)') "   atol = ", atol
    write (6, '(A,es9.2)') "   alpha = ", alpha
    write (6, *) "   ydot_i = -alpha*i * y_i (i = 1,...,neq)"
    write (6, *) "  "
  end if

  ! Create solution vector, point at its data, and set initial condition
  sunvec_y => FN_VNew_Parallel(comm, nlocal, neq, sunctx)
  y(1:nlocal) => FN_VGetArrayPointer(sunvec_y)
  y = 1.d0

  ! Create the ERKStep timestepper module
  arkode_mem = FERKStepCreate(c_funloc(frhs), t0, sunvec_y, sunctx)
  if (.not. c_associated(arkode_mem)) then
    print *, "Error: FERKStepCreate returned NULL"
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify tolerances
  retval = FARKodeSStolerances(arkode_mem, rtol, atol)
  if (retval /= 0) then
    print *, "Error: FARKodeSStolerances returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Main time-stepping loop: calls FARKodeEvolve to perform the
  ! integration, then prints results.  Stops when the final time
  ! has been reached.
  t(1) = T0
  dTout = 0.1d0
  tout = T0 + dTout
  if (outproc) then
    write (6, *) "        t      steps     steps att.    fe"
    write (6, *) "   -----------------------------------------"
  end if
  do ioutput = 1, Nt

    ! Integrate to output time
    retval = FARKodeEvolve(arkode_mem, tout, sunvec_y, t, ARK_NORMAL)
    if (retval /= 0) then
      print *, "Error: FARKodeEvolve returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FARKodeGetNumSteps(arkode_mem, nst)
    if (retval /= 0) then
      print *, "Error: FARKodeGetNumSteps returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
    if (retval /= 0) then
      print *, "Error: FARKodeGetNumStepAttempts returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FARKodeGetNumRhsEvals(arkode_mem, 0, nfe)
    if (retval /= 0) then
      print *, "Error: FARKodeGetNumRhsEvals returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    ! print solution stats and update internal time
    if (outproc) write (6, '(3x,f10.6,3(3x,i5))') t, nst, nst_a, nfe
    tout = min(tout + dTout, Tf)

  end do
  if (outproc) then
    write (6, *) "   -----------------------------------------"
  end if

  ! Get max. absolute error in the local vector.
  errmax = 0.d0
  do i = 1, nlocal
    erri = y(i) - exp(-alpha*(myid*nlocal + i)*t(1))
    errmax = max(errmax, abs(erri))
  end do

  ! Get global max. error from MPI_Reduce call.
  call MPI_Reduce(errmax, gerrmax, 1, MPI_DOUBLE, MPI_MAX, &
                  0, comm, ierr)
  if (ierr /= MPI_SUCCESS) then
    print *, "Error in MPI_Reduce = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Print global max. error
  if (outproc) print '(a,es10.2)', "Max. absolute error is ", gerrmax

  ! Get final statistics
  retval = FARKodeGetNumSteps(arkode_mem, nst)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumSteps returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumStepAttempts returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumRhsEvals(arkode_mem, 0, nfe)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumRhsEvals returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumErrTestFails(arkode_mem, netf)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumErrTestFails returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Print some final statistics
  if (outproc) then
    write (6, *) "  "
    write (6, *) "Final Solver Statistics:"
    write (6, '(2(A,i6),A)') "   Internal solver steps = ", nst, &
      " (attempted = ", nst_a, ")"
    write (6, '(A,i6)') "   Total RHS evals = ", nfe
    write (6, '(A,i6)') "   Total number of error test failures = ", netf
  end if

  ! Clean up and return with successful completion
  call FARKodeFree(arkode_mem)        ! free integrator memory
  call FN_VDestroy(sunvec_y)          ! free vector memory
  call MPI_Barrier(comm, ierr)
  call MPI_Finalize(ierr)             ! Finalize MPI

end program driver
!-----------------------------------------------------------------
