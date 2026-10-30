! ----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ UMBC
!                Radu Serban and Alan C. Hindmarsh @ LLNL
!                modified by Daniel M. Margolis @ UMBC
! ----------------------------------------------------------------
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
! ----------------------------------------------------------------
! Example problem for FIDA: 2D heat equation, parallel, GMRES,
! IDABBDPRE.
!
! This example solves a discretized 2D heat equation problem.
! This version uses the Krylov solver SUNSPGMR and BBD
! preconditioning.
!
! The DAE system solved is a spatial discretization of the PDE
!          du/dt = d^2u/dx^2 + d^2u/dy^2
! on the unit square. The boundary condition is u = 0 on all edges.
! Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
! PDE is treated with central differences on a uniform MX x MY
! grid. The values of u at the interior points satisfy ODEs, and
! equations u = 0 at the boundaries are appended, to form a DAE
! system of size N = MX * MY. Here MX = MY = 10.
!
! The system is actually implemented on submeshes, processor by
! processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
! processors.
!
! The system is solved with FIDA using the Krylov linear solver
! SUNSPGMR in conjunction with the preconditioner module IDABBDPRE.
! The preconditioner uses a tridiagonal approximation
! (half-bandwidths = 1). The constraints u >= 0 are posed for all
! components. Local error testing on the boundary values is
! suppressed. Output is taken at t = 0, .01, .02, .04, ..., 10.24.
! ----------------------------------------------------------------

module Heat2DKryBBD_mod
  ! --------------------------------------------------------------
  ! Description:
  !    Module containing problem-defining parameters, as well as
  !    data buffers for MPI exchanges with neighboring processes.
  !    Also contains routines to:
  !      (a) initialize the module
  !      (b) perform exchanges
  !      (c) free module data.
  ! --------------------------------------------------------------

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

  ! Grid and MPI domain decomposition information
  integer :: nx            ! global number of x grid points
  integer :: ny            ! global number of y grid points
  integer :: is            ! global x indices of this subdomain
  integer :: ie
  integer :: js            ! global y indices of this subdomain
  integer :: je
  integer :: nxl           ! local number of x grid points
  integer :: nyl           ! local number of y grid points
  integer(c_int64_t) :: N, Ntot
  real(c_double)  :: dx    ! x-directional mesh spacing
  real(c_double)  :: dy    ! y-directional mesh spacing
  integer, target :: comm  ! communicator object
  integer :: myid          ! MPI process ID
  integer :: nprocs        ! total number of MPI processes
  logical :: HaveNbor(2, 2) ! flags denoting neighbor on boundary
  real(c_double), dimension(:), allocatable :: Erecv  ! receive buffers for neighbor exchange
  real(c_double), dimension(:), allocatable :: Wrecv
  real(c_double), dimension(:), allocatable :: Nrecv
  real(c_double), dimension(:), allocatable :: Srecv
  real(c_double), dimension(:), allocatable :: Esend  ! send buffers for neighbor exchange
  real(c_double), dimension(:), allocatable :: Wsend
  real(c_double), dimension(:), allocatable :: Nsend
  real(c_double), dimension(:), allocatable :: Ssend

  ! Problem parameters
  integer(c_int64_t)  :: mudq, mldq, mu, ml
  integer(c_int)   :: maxl
  real(c_double) :: kx   ! x-directional diffusion coefficient
  real(c_double) :: ky   ! y-directional diffusion coefficient

  ! Printed parameters
  real(c_double)  :: h(1)       ! size of last step
  integer(c_int)  :: k(1)       ! number of last order
  integer(c_long) :: nst(1)     ! number of time steps
  integer(c_long) :: nre(1)     ! number of residual evals
  integer(c_long) :: nreLS(1)   ! number of residual linear solver evals
  integer(c_long) :: nge(1)     ! number of implicit RHS evals
  integer(c_long) :: netf(1)    ! number of error test fails
  integer(c_long) :: nncf(1)    ! number of nonlinear convergence fails
  integer(c_long) :: nlcf(1)    ! number of linear convergence fails
  integer(c_long) :: nni(1)     ! number of nonlinear iters
  integer(c_long) :: nli(1)     ! number of linear iters
  integer(c_long) :: npre(1)    ! number of preconditioner setups
  integer(c_long) :: npsol(1)   ! number of preconditioner solves

contains

  ! --------------------------------------------------------------
  ! Initialize memory allocated within Userdata (set to defaults)
  ! --------------------------------------------------------------
  subroutine InitHeat2DData()
    implicit none
    nx = 0
    ny = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    nxl = 0
    nyl = 0
    dx = 0.d0
    dy = 0.d0
    kx = 0.d0
    ky = 0.d0
    comm = MPI_COMM_WORLD
    myid = 0
    nprocs = 0
    HaveNbor = .false.
    if (allocated(Erecv)) deallocate (Erecv)
    if (allocated(Wrecv)) deallocate (Wrecv)
    if (allocated(Nrecv)) deallocate (Nrecv)
    if (allocated(Srecv)) deallocate (Srecv)
    if (allocated(Esend)) deallocate (Esend)
    if (allocated(Wsend)) deallocate (Wsend)
    if (allocated(Nsend)) deallocate (Nsend)
    if (allocated(Ssend)) deallocate (Ssend)
  end subroutine InitHeat2DData
  ! --------------------------------------------------------------

  ! --------------------------------------------------------------
  ! Set up parallel decomposition
  ! --------------------------------------------------------------
  subroutine SetupDecomp(ierr)
    ! declarations
    implicit none

    integer, intent(out) :: ierr
    integer :: dims(2), periods(2), coords(2)

    ! internals

    ! get suggested parallel decomposition
    dims = (/0, 0/)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in MPI_Comm_size = ", ierr
      return
    end if
    call MPI_Dims_create(nprocs, 2, dims, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in MPI_Dims_create = ", ierr
      return
    end if

    ! set up 2D Cartesian communicator
    periods = (/0, 0/)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in MPI_Cart_create = ", ierr
      return
    end if
    call MPI_Comm_rank(comm, myid, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in MPI_Comm_rank = ", ierr
      return
    end if

    ! determine local extents
    call MPI_Cart_get(comm, 2, dims, periods, coords, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in MPI_Cart_get = ", ierr
      return
    end if
    is = nx*coords(1)/dims(1) + 1
    ie = nx*(coords(1) + 1)/dims(1)
    js = ny*coords(2)/dims(2) + 1
    je = ny*(coords(2) + 1)/dims(2)
    nxl = ie - is + 1
    nyl = je - js + 1

    ! determine if I have neighbors, and allocate exchange buffers
    HaveNbor(1, 1) = (is /= 1)
    HaveNbor(1, 2) = (ie /= nx)
    HaveNbor(2, 1) = (js /= 1)
    HaveNbor(2, 2) = (je /= ny)
    if (HaveNbor(1, 1)) then
      allocate (Wrecv(nyl))
      allocate (Wsend(nyl))
    end if
    if (HaveNbor(1, 2)) then
      allocate (Erecv(nyl))
      allocate (Esend(nyl))
    end if
    if (HaveNbor(2, 1)) then
      allocate (Srecv(nxl))
      allocate (Ssend(nxl))
    end if
    if (HaveNbor(2, 2)) then
      allocate (Nrecv(nxl))
      allocate (Nsend(nxl))
    end if

    ierr = 0     ! return with success flag
    return
  end subroutine SetupDecomp
  ! --------------------------------------------------------------

  ! --------------------------------------------------------------
  ! Free memory allocated within Userdata
  ! --------------------------------------------------------------
  subroutine FreeHeat2DData(ierr)
    implicit none
    integer, intent(out) :: ierr
    if (allocated(Wrecv)) deallocate (Wrecv)
    if (allocated(Wsend)) deallocate (Wsend)
    if (allocated(Erecv)) deallocate (Erecv)
    if (allocated(Esend)) deallocate (Esend)
    if (allocated(Srecv)) deallocate (Srecv)
    if (allocated(Ssend)) deallocate (Ssend)
    if (allocated(Nrecv)) deallocate (Nrecv)
    if (allocated(Nsend)) deallocate (Nsend)
    ierr = 0     ! return with success flag
    return
  end subroutine FreeHeat2DData
  ! --------------------------------------------------------------

  subroutine InitProfile(sunvec_y, sunvec_ydot, sunvec_id, &
                         sunvec_res, sunvec_c, ierr)
    use fnvector_parallel_mod
    implicit none
    type(N_Vector), pointer, intent(inout) :: sunvec_y
    type(N_Vector), pointer, intent(inout) :: sunvec_ydot
    type(N_Vector), pointer, intent(inout) :: sunvec_id
    type(N_Vector), pointer, intent(inout) :: sunvec_res
    type(N_Vector), pointer, intent(inout) :: sunvec_c
    integer(c_int), intent(in) :: ierr
    real(c_double), pointer, dimension(nxl, nyl) :: y(:, :), ydot(:, :), id(:, :), res(:, :), cstr(:, :)
    real(c_double) :: xreal, yreal
    integer(c_int) :: retval
    type(c_ptr) :: user_data
    integer :: i, j

    user_data = c_null_ptr

    ! Create solution vector, point at its data, and set initial condition
    N = nxl*nyl
    Ntot = nx*ny
    sunvec_y => FN_VNew_Parallel(comm, N, Ntot, sunctx)
    sunvec_ydot => FN_VNew_Parallel(comm, N, Ntot, sunctx)
    sunvec_id => FN_VNew_Parallel(comm, N, Ntot, sunctx)
    sunvec_res => FN_VNew_Parallel(comm, N, Ntot, sunctx)
    sunvec_c => FN_VNew_Parallel(comm, N, Ntot, sunctx)
    y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)
    ydot(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_ydot)
    id(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_id)
    res(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_res)
    cstr(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_c)
    id = 1.d0
    do i = 1, nxl
      xreal = dx*dble(is + i - 2)
      do j = 1, nyl
        yreal = dy*dble(js + j - 2)
        if (.not. HaveNbor(1, 1) .and. i == 1) then
          id(i, j) = 0.d0
        end if
        if (.not. HaveNbor(1, 2) .and. i == nxl) then
          id(i, j) = 0.d0
        end if
        if (.not. HaveNbor(2, 1) .and. j == 1) then
          id(i, j) = 0.d0
        end if
        if (.not. HaveNbor(2, 2) .and. j == nyl) then
          id(i, j) = 0.d0
        end if
        y(i, j) = 16.d0*xreal*(1.d0 - xreal)*yreal*(1.d0 - yreal)
      end do
    end do
    ydot = 0.d0
    cstr = 1.d0
    retval = resfn(0.d0, sunvec_y, sunvec_ydot, sunvec_res, user_data)
    if (retval /= 0) then
      print *, "Error: resfn in InitProfile returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if
    ydot = -1.d0*res

    return
  end subroutine InitProfile
  ! --------------------------------------------------------------

  subroutine getStats(ida_mem, retval, ierr)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fida_mod
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr), intent(in)  :: ida_mem
    integer(c_int), intent(in)  :: ierr
    integer(c_int), intent(out) :: retval

    retval = FIDAGetLastOrder(ida_mem, k)
    if (retval /= 0) then
      print *, "Error: FIDAGetLastOrder returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumSteps(ida_mem, nst)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumSteps returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumNonlinSolvIters(ida_mem, nni)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumNonlinSolvIters returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumLinIters(ida_mem, nli)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumLinIters returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumResEvals(ida_mem, nre)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumResEvals returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumLinResEvals(ida_mem, nreLS)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumLinResEvals returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDABBDPrecGetNumGfnEvals(ida_mem, nge)
    if (retval /= 0) then
      print *, "Error: FIDABBDPrecGetNumGfnEvals returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetLastStep(ida_mem, h)
    if (retval /= 0) then
      print *, "Error: FIDAGetLastStep returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumPrecEvals(ida_mem, npre)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumPrecEvals returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumPrecSolves(ida_mem, npsol)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumPrecSolves returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

  end subroutine getStats
  ! --------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! DAE residual function r(t,y).
  ! ----------------------------------------------------------------
  integer(c_int) function resfn(t, sunvec_y, sunvec_ydot, sunvec_res, &
                                user_data) result(retval) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t            ! current time
    type(N_Vector)        :: sunvec_y     ! solution N_Vector
    type(N_Vector)        :: sunvec_ydot  ! rhs N_Vector
    type(N_Vector)        :: sunvec_res   ! residual N_Vector
    type(c_ptr), value :: user_data    ! user-defined data

    !======= Internals ============

    ! Exchange boundary data with neighbors
    retval = Exchange(N, t, sunvec_y, sunvec_ydot, sunvec_res, user_data)
    if (retval /= MPI_SUCCESS) then
      write (0, *) "Error in Exchange = ", retval
      return
    end if

    retval = LocalFn(N, t, sunvec_y, sunvec_ydot, sunvec_res, user_data)
    if (retval /= MPI_SUCCESS) then
      write (0, *) "Error in LocalFn = ", retval
      return
    end if

    retval = 0              ! Return with success
    return
  end function resfn
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Perform neighbor exchange (Communication function)
  ! ----------------------------------------------------------------
  integer(c_int) function Exchange(Nloc, t, sunvec_y, sunvec_ydot, &
                                   sunvec_g, user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    integer(c_int64_t), value :: Nloc
    real(c_double), value :: t            ! current time
    type(N_Vector)         :: sunvec_y     ! solution N_Vector
    type(N_Vector)         :: sunvec_ydot  ! rhs N_Vector
    type(N_Vector)         :: sunvec_g     ! evaluated N_Vector
    type(c_ptr), value :: user_data    ! user-defined data

    real(c_double), pointer, dimension(nxl, nyl) :: y(:, :)
    integer :: reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN; 
    integer :: stat(MPI_STATUS_SIZE)
    integer :: i, ipW, ipE, ipS, ipN
    integer :: coords(2), dims(2), periods(2), nbcoords(2)

    ! internals
    y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)

    ! MPI neighborhood information
    call MPI_Cart_get(comm, 2, dims, periods, coords, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in MPI_Cart_get = ", ierr
      return
    end if
    if (HaveNbor(1, 1)) then
      nbcoords = (/coords(1) - 1, coords(2)/)
      call MPI_Cart_rank(comm, nbcoords, ipW, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Cart_rank = ", ierr
        return
      end if
    end if
    if (HaveNbor(1, 2)) then
      nbcoords = (/coords(1) + 1, coords(2)/)
      call MPI_Cart_rank(comm, nbcoords, ipE, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Cart_rank = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 1)) then
      nbcoords = (/coords(1), coords(2) - 1/)
      call MPI_Cart_rank(comm, nbcoords, ipS, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Cart_rank = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 2)) then
      nbcoords = (/coords(1), coords(2) + 1/)
      call MPI_Cart_rank(comm, nbcoords, ipN, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Cart_rank = ", ierr
        return
      end if
    end if

    ! open Irecv buffers
    if (HaveNbor(1, 1)) then
      call MPI_Irecv(Wrecv, nyl, MPI_DOUBLE_PRECISION, ipW, &
                     MPI_ANY_TAG, comm, reqRW, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Irecv = ", ierr
        return
      end if
    end if
    if (HaveNbor(1, 2)) then
      call MPI_Irecv(Erecv, nyl, MPI_DOUBLE_PRECISION, ipE, &
                     MPI_ANY_TAG, comm, reqRE, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Irecv = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 1)) then
      call MPI_Irecv(Srecv, nxl, MPI_DOUBLE_PRECISION, ipS, &
                     MPI_ANY_TAG, comm, reqRS, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Irecv = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 2)) then
      call MPI_Irecv(Nrecv, nxl, MPI_DOUBLE_PRECISION, ipN, &
                     MPI_ANY_TAG, comm, reqRN, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Irecv = ", ierr
        return
      end if
    end if

    ! send data
    if (HaveNbor(1, 1)) then
      do i = 1, nyl
        Wsend(i) = y(1, i)
      end do
      call MPI_Isend(Wsend, nyl, MPI_DOUBLE_PRECISION, ipW, 0, &
                     comm, reqSW, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Isend = ", ierr
        return
      end if
    end if
    if (HaveNbor(1, 2)) then
      do i = 1, nyl
        Esend(i) = y(nxl, i)
      end do
      call MPI_Isend(Esend, nyl, MPI_DOUBLE_PRECISION, ipE, 1, &
                     comm, reqSE, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Isend = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 1)) then
      do i = 1, nxl
        Ssend(i) = y(i, 1)
      end do
      call MPI_Isend(Ssend, nxl, MPI_DOUBLE_PRECISION, ipS, 2, &
                     comm, reqSS, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Isend = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 2)) then
      do i = 1, nxl
        Nsend(i) = y(i, nyl)
      end do
      call MPI_Isend(Nsend, nxl, MPI_DOUBLE_PRECISION, ipN, 3, &
                     comm, reqSN, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Isend = ", ierr
        return
      end if
    end if

    ! wait for messages to finish
    if (HaveNbor(1, 1)) then
      call MPI_Wait(reqRW, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
      call MPI_Wait(reqSW, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
    end if
    if (HaveNbor(1, 2)) then
      call MPI_Wait(reqRE, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
      call MPI_Wait(reqSE, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 1)) then
      call MPI_Wait(reqRS, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
      call MPI_Wait(reqSS, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
    end if
    if (HaveNbor(2, 2)) then
      call MPI_Wait(reqRN, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
      call MPI_Wait(reqSN, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
        write (0, *) "Error in MPI_Wait = ", ierr
        return
      end if
    end if

    ierr = MPI_SUCCESS    ! return with success flag
    return
  end function Exchange
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Processor-local portion of the DAE residual function.
  ! ----------------------------------------------------------------
  integer(c_int) function LocalFn(Nloc, t, sunvec_y, sunvec_ydot, sunvec_g, &
                                  user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    integer(c_int64_t), value :: Nloc
    real(c_double), value :: t            ! current time
    type(N_Vector)         :: sunvec_y     ! solution N_Vector
    type(N_Vector)         :: sunvec_ydot  ! rhs N_Vector
    type(N_Vector)         :: sunvec_g     ! evaluated N_Vector
    type(c_ptr), value :: user_data    ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nxl, nyl) :: y(:, :)
    real(c_double), pointer, dimension(nxl, nyl) :: ydot(:, :)
    real(c_double), pointer, dimension(nxl, nyl) :: res(:, :)

    ! local data
    real(c_double) :: c1, c2, c3
    integer :: i, j

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)
    ydot(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_ydot)
    res(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_g)

    ! set constants
    c1 = kx/dx/dx
    c2 = ky/dy/dy
    c3 = -2.d0*(c1 + c2)

    ! Copy y NVector into res NVector
    res = y

    ! iterate over subdomain boundaries (if not at overall domain boundary)
    do i = 1, nxl
      do j = 1, nyl
        if (i == 1 .and. j == 1) then
          if (HaveNbor(1, 1) .and. HaveNbor(2, 1)) then  ! South-West corner
            res(i, j) = c1*(Wrecv(j) + y(i + 1, j)) + c2*(Srecv(i) + y(i, j + 1)) + c3*y(i, j)
          end if
        else if (i == 1 .and. j == nyl) then
          if (HaveNbor(1, 1) .and. HaveNbor(2, 2)) then  ! North-West corner
            res(i, j) = c1*(Wrecv(j) + y(i + 1, j)) + c2*(y(i, j - 1) + Nrecv(i)) + c3*y(i, j)
          end if
        else if (i == nxl .and. j == 1) then
          if (HaveNbor(1, 2) .and. HaveNbor(2, 1)) then  ! South-East corner
            res(i, j) = c1*(y(i - 1, j) + Erecv(j)) + c2*(Srecv(i) + y(i, j + 1)) + c3*y(i, j)
          end if
        else if (i == nxl .and. j == nyl) then
          if (HaveNbor(1, 2) .and. HaveNbor(2, 2)) then  ! North-East corner
            res(i, j) = c1*(y(i - 1, j) + Erecv(j)) + c2*(y(i, j - 1) + Nrecv(i)) + c3*y(i, j)
          end if
        else if (i == 1) then
          if (HaveNbor(1, 1)) then                      ! West face
            res(i, j) = c1*(Wrecv(j) + y(i + 1, j)) + c2*(y(i, j - 1) + y(i, j + 1)) + c3*y(i, j)
          end if
        else if (i == nxl) then
          if (HaveNbor(1, 2)) then                      ! East face
            res(i, j) = c1*(y(i - 1, j) + Erecv(j)) + c2*(y(i, j - 1) + y(i, j + 1)) + c3*y(i, j)
          end if
        else if (j == 1) then
          if (HaveNbor(2, 1)) then                      ! South face
            res(i, j) = c1*(y(i - 1, j) + y(i + 1, j)) + c2*(Srecv(i) + y(i, j + 1)) + c3*y(i, j)
          end if
        else if (j == nyl) then
          if (HaveNbor(2, 2)) then                      ! North face
            res(i, j) = c1*(y(i - 1, j) + y(i + 1, j)) + c2*(y(i, j - 1) + Nrecv(i)) + c3*y(i, j)
          end if
        else
          res(i, j) = c1*(y(i - 1, j) + y(i + 1, j)) + c2*(y(i, j - 1) + y(i, j + 1)) + c3*y(i, j)
        end if
        res(i, j) = ydot(i, j) - res(i, j)
      end do
    end do

    ierr = 0     ! Return with success
    return
  end function LocalFn
  ! ----------------------------------------------------------------

end module Heat2DKryBBD_mod
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program driver

  ! inclusions
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fida_mod                   ! Access IDA
  use fnvector_parallel_mod      ! Access parallel N_Vector
  use fsunlinsol_spgmr_mod       ! Access GMRES SUNLinearSolver
  use Heat2DKryBBD_mod

  !======= Declarations =========
  implicit none

  ! Declarations
  ! general problem parameters
  integer, parameter :: Nt = 11           ! total number of output times
  integer, parameter :: nx_ = 100         ! spatial mesh size
  integer, parameter :: ny_ = 100
  real(c_double), parameter :: T0 = 0.d0        ! initial time
  real(c_double), parameter :: T1 = 0.01d0      ! final time
  real(c_double), parameter :: rtol = 0.d0      ! relative and absolute tolerances
  real(c_double), parameter :: atol = 1.d-3
  real(c_double), parameter :: kx_ = 1.0d0      ! heat conductivity coefficients
  real(c_double), parameter :: ky_ = 1.0d0

  ! solution vector and other local variables
  type(N_Vector), pointer :: sunvec_y        ! solution N_Vector
  type(N_Vector), pointer :: sunvec_f        ! derivative N_Vector
  type(N_Vector), pointer :: sunvec_id       ! derivative N_Vector
  type(N_Vector), pointer :: sunvec_res      ! derivative N_Vector
  type(N_Vector), pointer :: sunvec_c        ! constraint N_Vector
  real(c_double), pointer, dimension(nxl, nyl) :: y(:, :)  ! vector data
  type(SUNLinearSolver), pointer :: sun_LS   ! linear solver
  type(SUNMatrix), pointer :: sunmat_A ! sundials matrix
  type(c_ptr)     :: ida_mem                 ! IDA memory
  integer(c_int) :: retval
  integer :: ierr, case
  logical :: outproc
  real(c_double) :: t(1), tout, ymax
  integer :: i, j, ioutput
  character(100) :: outname

  ! initialize MPI
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Init = ", ierr
    stop 1
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Comm_rank = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Initialize Heat2DData module
  call InitHeat2DData()
  nx = nx_
  ny = ny_
  kx = kx_
  ky = ky_
  dx = 1.d0/dble(nx - 1)   ! x mesh spacing
  dy = 1.d0/dble(ny - 1)   ! x mesh spacing

  ! Set up parallel decomposition (computes local mesh sizes)
  call SetupDecomp(ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in SetupDecomp = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if

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
    write (6, *) "2D Heat PDE test problem:"; 
    write (6, '(A,i4)') "   nprocs = ", nprocs
    write (6, '(A,i4)') "   nx = ", nx
    write (6, '(A,i4)') "   ny = ", ny
    write (6, '(A,f5.2)') "   kx = ", kx
    write (6, '(A,f5.2)') "   ky = ", ky
    write (6, '(A,es9.2)') "   rtol = ", rtol
    write (6, '(A,es9.2)') "   atol = ", atol
    write (6, '(A,i4)') "   nxl (proc 0) = ", nxl
    write (6, '(A,i4)') "   nyl (proc 0) = ", nyl
    write (6, *) "  "
  end if

  ! Create the IDA timestepper module
  ida_mem = FIDACreate(sunctx)
  if (.not. c_associated(ida_mem)) then
    print *, "Error: FIDACreate returned NULL"
    call MPI_Abort(comm, 1, ierr)
  end if

  call InitProfile(sunvec_y, sunvec_f, sunvec_id, sunvec_res, sunvec_c, ierr)

  retval = FIDAInit(ida_mem, c_funloc(resfn), t0, sunvec_y, sunvec_f)
  if (retval /= 0) then
    print *, "Error: FIDAInit returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Create linear solver
  maxl = 0
  sun_LS => FSUNLinSol_SPGMR(sunvec_y, SUN_PREC_LEFT, maxl, sunctx)
  if (.not. associated(sun_LS)) then
    print *, "Error: FSUNLinSol_PCG returned NULL"
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Attach linear solver
  sunmat_A => null()
  retval = FIDASetLinearSolver(ida_mem, sun_LS, sunmat_A)
  if (retval /= 0) then
    print *, "Error: FIDASetLinearSolver returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FSUNLinSol_SPGMRSetMaxRestarts(sun_LS, 5)
  if (retval /= 0) then
    print *, "Error: FSUNLinSol_SPGMRSetMaxRestarts returned", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Attach preconditioner
  mudq = nxl
  mldq = nxl
  mu = 1
  ml = 1
  retval = FIDABBDPrecInit(ida_mem, N, mudq, mldq, mu, ml, 0.d0, &
                           c_funloc(LocalFn), c_funloc(Exchange))
  if (retval /= 0) then
    print *, "Error: FIDASetPreconditioner returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify tolerances
  retval = FIDASStolerances(ida_mem, rtol, atol)
  if (retval /= 0) then
    print *, "Error: FIDASStolerances returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FIDASetSuppressAlg(ida_mem, SUNTRUE)
  if (retval /= 0) then
    print *, "Error: FIDASetSuppressAlg returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FIDASetId(ida_mem, sunvec_id)
  if (retval /= 0) then
    print *, "Error: FIDASetId returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FIDASetConstraints(ida_mem, sunvec_c)
  if (retval /= 0) then
    print *, "Error: FIDASetConstraints returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Each processor outputs subdomain information
  write (outname, '(16Hheat2d_subdomain,f4.3,4H.txt)') myid/1000.0
  open (100, file=outname)
  write (100, '(6(i9,1x))') nx, ny, is, ie, js, je
  close (100)

  ! Open output streams for results, access data array
  write (outname, '(6Hheat2d,f4.3,4H.txt)') myid/1000.0
  open (101, file=outname)

  ! Output initial condition to disk
  y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)
  do j = 1, nyl
    do i = 1, nxl
      write (101, '(es25.16)', advance='no') y(i, j)
    end do
  end do
  write (101, *) "  "

  do case = 1, 2
    if (case == 2) then
      mudq = 5
      mldq = 5
      call InitProfile(sunvec_y, sunvec_f, sunvec_id, sunvec_res, sunvec_c, ierr)
      retval = FIDAReInit(ida_mem, t0, sunvec_y, sunvec_f)
      if (retval /= 0) then
        print *, "Error: FIDAReInit returned ", retval
        call MPI_Abort(comm, 1, ierr)
      end if
      retval = FIDABBDPrecReInit(ida_mem, mudq, mldq, 0.d0)
      if (retval /= 0) then
        print *, "Error: FIDABBDPrecReInit returned ", retval
        call MPI_Abort(comm, 1, ierr)
      end if
      if (outproc) then
        write (6, *) "  "
        write (6, *) "Case ", case
        write (6, *) "  Difference quotient half-bandwidths = ", mudq
        write (6, *) "  Retained matrix half-bandwidths = ", mu
        write (6, *) "  "
        write (6, *) "Output Summary"
        write (6, *) "  umax = max-norm of solution"
        write (6, *) "  nre = nre + nreLS (total number of RES evals.)"
      end if
    end if
    if (case == 1) then
      if (outproc) then
        write (6, *) "  "
        write (6, *) "Case ", case
        write (6, *) "  Difference quotient half-bandwidths = ", mudq
        write (6, *) "  Retained matrix half-bandwidths = ", mu
        write (6, *) "  "
        write (6, *) "Output Summary"
        write (6, *) "  umax = max-norm of solution"
        write (6, *) "  nre = nre + nreLS (total number of RES evals.)"
      end if
    end if
    ! Main time-stepping loop: calls IDA to perform the integration, then
    ! prints results.  Stops when the final time has been reached
    t = T0
    tout = T1
    if (outproc) then
      write (6, *) "   "
      write (6, *) "        t      ||u||_max       k  nst   nni   nli   nre    nge       h      npe   nps"
      write (6, *) "   ------------------------------------------------------------------------------------"
    end if
    do ioutput = 1, Nt

      ! Integrate to output time
      retval = FIDASolve(ida_mem, tout, t, sunvec_y, sunvec_f, IDA_NORMAL)
      if (retval /= 0) then
        print *, "Error: FIDASolve returned ", retval
        call MPI_Abort(comm, 1, ierr)
      end if

      ! print solution stats and update internal time
      ymax = FN_VMaxNorm(sunvec_y)
      call getStats(ida_mem, retval, ierr)
      if (outproc) write (6, '(2x,f10.6,2x,es13.5,3x,i1,3x,i2,3x,i3,3x,i3,3x,i3,a,i3,3x,i4,3x,es9.2,3x,i2,3x,i3)') &
        t, ymax, k, nst, nni, nli, nre, "+", nreLS, nge, h, npre, npsol
      tout = 2.0d0*tout

      ! output results to disk
      do j = 1, nyl
        do i = 1, nxl
          write (101, '(es25.16)', advance='no') y(i, j)
        end do
      end do
      write (101, *) "  "

    end do
    if (outproc) then
      write (6, *) "   ------------------------------------------------------------------------------------"
    end if
    close (101)

    retval = FIDAGetNumErrTestFails(ida_mem, netf)
    if (retval /= 0) then
      print *, "Error: FIDAGetNumErrTestFails returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumNonlinSolvConvFails(ida_mem, nncf)
    if (retval /= 0) then
      print *, "Error: FIDAInit returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    retval = FIDAGetNumLinConvFails(ida_mem, nlcf)
    if (retval /= 0) then
      print *, "Error: FIDAInit returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    ! Print some final statistics
    if (outproc) then
      write (6, *) "  "
      write (6, *) "Final Solver Statistics:"
      write (6, '(A,i6)') "   Total number of error test failures      = ", netf
      write (6, '(A,i6)') "   Total number of nonlinear conv. failures = ", nncf
      write (6, '(A,i6)') "   Total number of linear conv. failures    = ", nlcf
    end if
  end do

  ! Clean up and return with successful completion
  call FIDAFree(ida_mem)           ! free integrator memory
  call FN_VDestroy_Parallel(sunvec_y)          ! free vector memory
  call FN_VDestroy_Parallel(sunvec_f)
  call FN_VDestroy_Parallel(sunvec_id)
  call FN_VDestroy_Parallel(sunvec_res)
  call FN_VDestroy_Parallel(sunvec_c)
  retval = FSUNLinSolFree(sun_LS)     ! free linear solver
  call FreeHeat2DData(ierr)           ! free user data
  call MPI_Barrier(comm, ierr)
  call MPI_Finalize(ierr)             ! Finalize MPI

end program driver
! ----------------------------------------------------------------
