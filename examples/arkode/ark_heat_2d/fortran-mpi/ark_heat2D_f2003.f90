! ----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ UMBC
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
! Example problem:
!
! The following test simulates a simple anisotropic 2D heat
! equation,
!    u_t = kx*u_xx + ky*u_yy + h,
! for t in [0, 10], (x,y) in [0, 1]^2, with initial conditions
!    u(0,x,y) =  0,
! stationary boundary conditions, i.e.
!    u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
! and a heat source of the form
!    h(x,y) = sin(pi*x)*sin(2*pi*y).
!
! Under this setup, the problem has an analytical solution:
!    u(t,x,y) = a(t)*sin(pi*x)*sin(2*pi*y), where
!    a(t) = (1 - exp(-(kx+4*ky)*pi^2*t)) / ((kx+4*ky)*pi^2).
!
! The spatial derivatives are computed using second-order
! centered differences, with the data distributed over nx*ny
! points on a uniform spatial grid.
!
! The spatial grid parameters nx and ny, the parameters kx and ky,
! as well as the desired relative and absolute solver tolerances,
! are provided in the input file input_heat2D.txt.
!
! This program solves the problem with a DIRK method.  This
! employs a Newton iteration with the PCG iterative linear solver,
! which itself uses a Jacobi preconditioner.  The example uses the
! built-in finite-difference Jacobian-vector product routine, but
! supplies both the RHS and preconditioner setup/solve functions.
!
! 20 outputs are printed at equal intervals, and run statistics
! are printed at the end.
! ----------------------------------------------------------------

module Heat2DData
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
  real(c_double) :: dx     ! x-directional mesh spacing
  real(c_double) :: dy     ! y-directional mesh spacing
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
  real(c_double) :: kx   ! x-directional diffusion coefficient
  real(c_double) :: ky   ! y-directional diffusion coefficient
  real(c_double), dimension(:, :), allocatable :: h   ! heat source vector

  ! Preconditioning data
  real(c_double), dimension(:, :), allocatable :: d   ! inverse of Jacobian diagonal

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
    if (allocated(h)) deallocate (h)
    if (allocated(d)) deallocate (d)
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

    ! check that this has not been called before
    if (allocated(h) .or. allocated(d)) then
      write (0, *) "SetupDecomp warning: parallel decomposition already set up"
      ierr = 1
      return
    end if

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

    ! allocate temporary vectors
    allocate (h(nxl, nyl))    ! Create vector for heat source
    allocate (d(nxl, nyl))    ! Create vector for Jacobian diagonal

    ierr = 0     ! return with success flag
    return
  end subroutine SetupDecomp
  ! --------------------------------------------------------------

  ! --------------------------------------------------------------
  ! Perform neighbor exchange
  ! --------------------------------------------------------------
  subroutine Exchange(y, ierr)
    ! declarations
    implicit none
    real(c_double), intent(in) :: y(nxl, nyl)
    integer, intent(out) :: ierr
    integer :: reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN; 
    integer :: stat(MPI_STATUS_SIZE)
    integer :: i, ipW, ipE, ipS, ipN
    integer :: coords(2), dims(2), periods(2), nbcoords(2)

    ! internals

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
  end subroutine Exchange
  ! --------------------------------------------------------------

  ! --------------------------------------------------------------
  ! Free memory allocated within Userdata
  ! --------------------------------------------------------------
  subroutine FreeHeat2DData(ierr)
    implicit none
    integer, intent(out) :: ierr
    if (allocated(h)) deallocate (h)
    if (allocated(d)) deallocate (d)
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

  ! ----------------------------------------------------------------
  ! ODE RHS function f(t,y).
  ! ----------------------------------------------------------------
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
    real(c_double), pointer, dimension(nxl, nyl) :: y(:, :)
    real(c_double), pointer, dimension(nxl, nyl) :: ydot(:, :)

    ! local data
    real(c_double) :: c1, c2, c3
    integer :: i, j, ierr

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)
    ydot(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_ydot)

    ! Initialize ydot to zero
    ydot = 0.d0

    ! Exchange boundary data with neighbors
    call Exchange(y, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (0, *) "Error in Exchange = ", ierr
      retval = -1
      return
    end if

    ! iterate over subdomain interior, computing approximation to RHS
    c1 = kx/dx/dx
    c2 = ky/dy/dy
    c3 = -2.d0*(c1 + c2)
    do j = 2, nyl - 1
      do i = 2, nxl - 1
        ydot(i, j) = c1*(y(i - 1, j) + y(i + 1, j)) + c2*(y(i, j - 1) + y(i, j + 1)) + c3*y(i, j)
      end do
    end do

    ! iterate over subdomain boundaries (if not at overall domain boundary)
    if (HaveNbor(1, 1)) then    ! West face
      i = 1
      do j = 2, nyl - 1
        ydot(i, j) = c1*(Wrecv(j) + y(i + 1, j)) + c2*(y(i, j - 1) + y(i, j + 1)) + c3*y(i, j)
      end do
    end if
    if (HaveNbor(1, 2)) then    ! East face
      i = nxl
      do j = 2, nyl - 1
        ydot(i, j) = c1*(y(i - 1, j) + Erecv(j)) + c2*(y(i, j - 1) + y(i, j + 1)) + c3*y(i, j)
      end do
    end if
    if (HaveNbor(2, 1)) then    ! South face
      j = 1
      do i = 2, nxl - 1
        ydot(i, j) = c1*(y(i - 1, j) + y(i + 1, j)) + c2*(Srecv(i) + y(i, j + 1)) + c3*y(i, j)
      end do
    end if
    if (HaveNbor(2, 2)) then    ! West face
      j = nyl
      do i = 2, nxl - 1
        ydot(i, j) = c1*(y(i - 1, j) + y(i + 1, j)) + c2*(y(i, j - 1) + Nrecv(i)) + c3*y(i, j)
      end do
    end if
    if (HaveNbor(1, 1) .and. HaveNbor(2, 1)) then  ! South-West corner
      i = 1
      j = 1
      ydot(i, j) = c1*(Wrecv(j) + y(i + 1, j)) + c2*(Srecv(i) + y(i, j + 1)) + c3*y(i, j)
    end if
    if (HaveNbor(1, 1) .and. HaveNbor(2, 2)) then  ! North-West corner
      i = 1
      j = nyl
      ydot(i, j) = c1*(Wrecv(j) + y(i + 1, j)) + c2*(y(i, j - 1) + Nrecv(i)) + c3*y(i, j)
    end if
    if (HaveNbor(1, 2) .and. HaveNbor(2, 1)) then  ! South-East corner
      i = nxl
      j = 1
      ydot(i, j) = c1*(y(i - 1, j) + Erecv(j)) + c2*(Srecv(i) + y(i, j + 1)) + c3*y(i, j)
    end if
    if (HaveNbor(1, 2) .and. HaveNbor(2, 2)) then  ! North-East corner
      i = nxl
      j = nyl
      ydot(i, j) = c1*(y(i - 1, j) + Erecv(j)) + c2*(y(i, j - 1) + Nrecv(i)) + c3*y(i, j)
    end if

    ydot = ydot + h         ! add in heat source

    retval = 0              ! Return with success
    return
  end function frhs
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Preconditioner setup routine (fills inverse of Jacobian diagonal)
  ! ----------------------------------------------------------------
  integer(c_int) function PSetup(t, sunvec_y, sunvec_ydot, jok, jcurPtr, &
                                 gamma, user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t            ! current time
    type(N_Vector)        :: sunvec_y     ! solution N_Vector
    type(N_Vector)        :: sunvec_ydot  ! rhs N_Vector
    integer(c_int), value :: jok          ! flag to signal for Jacobian update
    integer(c_int)        :: jcurPtr      ! flag to signal Jacobian is current
    real(c_double), value :: gamma        ! current gamma value
    type(c_ptr), value :: user_data    ! user-defined data

    ! local variables
    real(c_double) :: c

    !======= Internals ============

    ! set constant for matrix diagonal
    c = 1.d0 + gamma*2.d0*(kx/dx/dx + ky/dy/dy)

    ! set all entries of d to the inverse of the diagonal values in interior
    ! (since boundary RHS is 0, set boundary diagonals to the same)
    d = 1.d0/c

    jcurPtr = 1  ! update jcur flag

    ierr = 0     ! Return with success
    return
  end function PSetup
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Preconditioner solve routine
  ! ----------------------------------------------------------------
  integer(c_int) function PSolve(t, sunvec_y, sunvec_ydot, sunvec_r, &
                                 sunvec_z, gamma, delta, lr, user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t            ! current time
    type(N_Vector)        :: sunvec_y     ! solution N_Vector
    type(N_Vector)        :: sunvec_ydot  ! rhs N_Vector
    type(N_Vector)        :: sunvec_r     ! rhs N_Vector
    type(N_Vector)        :: sunvec_z     ! rhs N_Vector
    real(c_double), value :: gamma        ! current gamma value
    real(c_double), value :: delta        ! current delta value
    integer(c_int), value :: lr           ! left or right preconditioning
    type(c_ptr), value :: user_data    ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nxl, nyl) :: r(:, :)
    real(c_double), pointer, dimension(nxl, nyl) :: z(:, :)

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    r(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_r)
    z(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_z)

    ! perform Jacobi solve (whole array operation)
    z = r*d

    ierr = 0     ! Return with success
    return
  end function PSolve
  ! ----------------------------------------------------------------

end module Heat2DData
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program driver

  ! inclusions
  use, intrinsic :: iso_c_binding
  use farkode_mod                ! Access ARKODE
  use farkode_arkstep_mod        ! Access ARKStep
  use fnvector_parallel_mod      ! Access parallel N_Vector
  use fsunlinsol_pcg_mod         ! Access GMRES SUNLinearSolver
  use Heat2DData

  !======= Declarations =========
  implicit none

  ! Declarations
  ! general problem parameters
  real(c_double), parameter :: pi = 3.1415926535897932d0
  integer, parameter :: Nt = 20           ! total number of output times
  integer, parameter :: nx_ = 60          ! spatial mesh size
  integer, parameter :: ny_ = 120
  real(c_double), parameter :: T0 = 0.d0        ! initial time
  real(c_double), parameter :: Tf = 0.3d0       ! final time
  real(c_double), parameter :: rtol = 1.d-5     ! relative and absolute tolerances
  real(c_double), parameter :: atol = 1.d-10
  real(c_double), parameter :: kx_ = 0.5d0      ! heat conductivity coefficients
  real(c_double), parameter :: ky_ = 0.75d0
  real(c_double), parameter :: nlscoef = 1.d-7  ! nonlinear solver tolerance factor

  ! solution vector and other local variables
  type(N_Vector), pointer :: sunvec_y                   ! solution N_Vector
  type(N_Vector), pointer :: sunvec_ones                ! masking vector for output
  real(c_double), pointer, dimension(nxl, nyl) :: y(:, :) ! vector data
  type(SUNLinearSolver), pointer :: sun_LS              ! linear solver
  type(SUNMatrix), pointer :: sunmat_A            ! sundials matrix
  type(c_ptr)     :: arkode_mem                         ! ARKODE memory
  integer(c_int64_t) :: N, Ntot
  integer(c_int) :: retval
  integer :: ierr
  logical :: outproc
  real(c_double) :: t(1), dTout, tout, urms
  integer(c_long) :: nst(1)     ! number of time steps
  integer(c_long) :: nst_a(1)   ! number of step attempts
  integer(c_long) :: nfe(1)     ! number of explicit RHS evals
  integer(c_long) :: nfi(1)     ! number of implicit RHS evals
  integer(c_long) :: netf(1)    ! number of error test fails
  integer(c_long) :: nni(1)     ! number of nonlinear iters
  integer(c_long) :: ncfn(1)    ! number of convergence fails
  integer(c_long) :: nli(1)     ! number of linear iters
  integer(c_long) :: npre(1)    ! number of preconditioner setups
  integer(c_long) :: npsol(1)   ! number of preconditioner solves
  integer :: i, j, ioutput
  character(len=4) :: idstring
  character(len=14) :: outname
  character(len=24) :: subdomainname

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
  dx = 1.d0/(nx - 1)   ! x mesh spacing
  dy = 1.d0/(ny - 1)   ! x mesh spacing

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

  ! Create solution vector, point at its data, and set initial condition
  N = nxl*nyl
  Ntot = nx*ny
  sunvec_y => FN_VNew_Parallel(comm, N, Ntot, sunctx)
  y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)
  y = 0.d0

  ! Create the ARKStep timestepper module
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(frhs), t0, sunvec_y, sunctx)
  if (.not. c_associated(arkode_mem)) then
    print *, "Error: FARKStepCreate returned NULL"
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Create linear solver
  sun_LS => FSUNLinSol_PCG(sunvec_y, SUN_PREC_LEFT, int(20, c_int), sunctx)
  if (.not. associated(sun_LS)) then
    print *, "Error: FSUNLinSol_PCG returned NULL"
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Attach linear solver
  sunmat_A => null()
  retval = FARKodeSetLinearSolver(arkode_mem, sun_LS, sunmat_A)
  if (retval /= 0) then
    print *, "Error: FARKodeSetLinearSolver returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Attach preconditioner
  retval = FARKodeSetPreconditioner(arkode_mem, c_funloc(PSetup), &
                                    c_funloc(PSolve))
  if (retval /= 0) then
    print *, "Error: FARKodeSetPreconditioner returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify tolerances
  retval = FARKodeSStolerances(arkode_mem, rtol, atol)
  if (retval /= 0) then
    print *, "Error: FARKodeSStolerances returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify nonlinear solver convergence coefficient
  retval = FARKodeSetNonlinConvCoef(arkode_mem, nlscoef)
  if (retval /= 0) then
    print *, "Error: FARKodeSetNonlinConvCoef returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify nonlinear solver predictor method
  retval = FARKodeSetPredictorMethod(arkode_mem, int(1, c_int))
  if (retval /= 0) then
    print *, "Error: FARKodeSetNonlinConvCoef returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify that problem is linearly implicit
  retval = FARKodeSetLinear(arkode_mem, int(0, c_int))
  if (retval /= 0) then
    print *, "Error: FARKodeSetNonlinConvCoef returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! fill in the heat source array
  do j = 1, nyl
    do i = 1, nxl
      h(i, j) = sin(pi*(is + i - 2)*dx)*sin(2.d0*pi*(js + j - 2)*dy)
    end do
  end do

  ! Each processor outputs subdomain information
  write (idstring, "(f4.3)") myid/1000.0
  subdomainname = "heat2d_subdomain"//idstring//".txt"
  open (100, file=subdomainname)
  write (100, '(6(i9,1x))') nx, ny, is, ie, js, je
  close (100)

  ! Open output streams for results, access data array
  outname = "heat2d"//idstring//".txt"
  open (101, file=outname)

  ! Output initial condition to disk
  do j = 1, nyl
    do i = 1, nxl
      write (101, '(es25.16)', advance='no') y(i, j)
    end do
  end do
  write (101, *) "  "

  ! create masking vector for computing solution RMS norm
  sunvec_ones => FN_VNew_Parallel(comm, N, Ntot, sunctx)
  call FN_VConst(1.d0, sunvec_ones)

  ! Main time-stepping loop: calls ARKODE to perform the integration, then
  ! prints results.  Stops when the final time has been reached
  t(1) = T0
  dTout = (Tf - T0)/Nt
  tout = T0 + dTout
  urms = FN_VWrmsNorm(sunvec_y, sunvec_ones)
  if (outproc) then
    write (6, *) "        t      ||u||_rms"
    write (6, *) "   ----------------------"
    write (6, '(2(2x,f10.6))') t, urms
  end if
  do ioutput = 1, Nt

    ! Integrate to output time
    retval = FARKodeEvolve(arkode_mem, tout, sunvec_y, t, ARK_NORMAL)
    if (retval /= 0) then
      print *, "Error: FARKodeEvolve returned ", retval
      call MPI_Abort(comm, 1, ierr)
    end if

    ! print solution stats and update internal time
    urms = FN_VWrmsNorm(sunvec_y, sunvec_ones)
    if (outproc) write (6, '(2(2x,f10.6))') t, urms
    tout = min(tout + dTout, Tf)

    ! output results to disk
    do j = 1, nyl
      do i = 1, nxl
        write (101, '(es25.16)', advance='no') y(i, j)
      end do
    end do
    write (101, *) "  "

  end do
  if (outproc) then
    write (6, *) "   ----------------------"
  end if
  close (101)

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

  retval = FARKodeGetNumRhsEvals(arkode_mem, 1, nfi)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumRhsEvals returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumErrTestFails(arkode_mem, netf)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumErrTestFails returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumNonlinSolvIters(arkode_mem, nni)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumNonlinSolvIters returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumLinConvFails(arkode_mem, ncfn)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumLinConvFails returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumLinIters(arkode_mem, nli)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumLinIters returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumPrecEvals(arkode_mem, npre)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumPrecEvals returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKodeGetNumPrecSolves(arkode_mem, npsol)
  if (retval /= 0) then
    print *, "Error: FARKodeGetNumPrecSolves returned ", retval
    call MPI_Abort(comm, 1, ierr)
  end if

  ! Print some final statistics
  if (outproc) then
    write (6, *) "  "
    write (6, *) "Final Solver Statistics:"
    write (6, '(2(A,i6),A)') "   Internal solver steps = ", nst, &
      " (attempted = ", nst_a, ")"
    write (6, '(A,i6,A,i6)') "   Total RHS evals:  Fe = ", nfe, ",  Fi = ", nfi
    write (6, '(A,i6)') "   Total linear iterations = ", nli
    write (6, '(A,i6)') "   Total number of Preconditioner setups = ", npre
    write (6, '(A,i6)') "   Total number of Preconditioner solves = ", npsol
    write (6, '(A,i6)') "   Total number of linear solver convergence failures = ", &
      ncfn
    write (6, '(A,i6)') "   Total number of Newton iterations = ", nni
    write (6, '(A,i6)') "   Total number of error test failures = ", netf
  end if

  ! Clean up and return with successful completion
  call FARKodeFree(arkode_mem)        ! free integrator memory
  call FN_VDestroy(sunvec_y)          ! free vector memory
  call FN_VDestroy(sunvec_ones)
  retval = FSUNLinSolFree(sun_LS)     ! free linear solver
  call FreeHeat2DData(ierr)           ! free user data
  call MPI_Barrier(comm, ierr)
  call MPI_Finalize(ierr)             ! Finalize MPI

end program driver
! ----------------------------------------------------------------
