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
! band SUNMatrix implementation.
! -----------------------------------------------------------------

module test_fsunmatrix_band
  use, intrinsic :: iso_c_binding
  use test_utilities
  implicit none

  integer(kind=myindextype), parameter :: N = 10
  integer(kind=myindextype), parameter :: mu = 2
  integer(kind=myindextype), parameter :: ml = 2

contains

  integer(c_int) function smoke_tests() result(fails)

    !======== Inclusions ==========
    use, intrinsic :: iso_c_binding

    use fsunmatrix_band_mod
    use fnvector_serial_mod

    !======== Declarations ========
    implicit none

    ! local variables
    type(SUNMatrix), pointer :: A, B               ! SUNMatrix
    type(N_Vector), pointer :: x, y               ! NVectors
    real(c_double), pointer :: matdat(:)          ! matrix data pointer
    integer(c_long)          :: lenrw(1), leniw(1) ! matrix real and int work space size
    integer(c_long)          :: val
    type(c_ptr), pointer :: cptr

    fails = 0
    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)

    !===== Calls to interface =====

    ! constructor
    A => FSUNBandMatrix(N, mu, ml, sunctx)
    if (.not. associated(A)) then
      print *, '>>> FAILED - ERROR in FSUNBandMatrix; halting'
      fails = 1
      return
    end if

    ! misc. matrix functions
    val = FSUNMatGetID_Band(A)
    val = FSUNBandMatrix_Rows(A)
    val = FSUNBandMatrix_Columns(A)
    val = FSUNBandMatrix_LowerBandwidth(A)
    val = FSUNBandMatrix_UpperBandwidth(A)
    val = FSUNBandMatrix_StoredUpperBandwidth(A)
    val = FSUNBandMatrix_LDim(A)
    matdat => FSUNBandMatrix_Data(A)
    cptr => FSUNBandMatrix_Cols(A)
    matdat => FSUNBandMatrix_Column(A, N)

    ! matrix operations
    B => FSUNMatClone_Band(A)
    if (.not. associated(B)) then
      print *, '>>> FAILED - ERROR in FSUNMatClone_Band; halting'
      fails = 1
      return
    end if
    fails = fails + FSUNMatZero_Band(A)
    fails = fails + FSUNMatCopy_Band(A, B)
    fails = fails + FSUNMatScaleAdd_Band(ONE, A, B)
    fails = fails + FSUNMatScaleAddI_Band(ONE, A)
    fails = fails + FSUNMatMatvec_Band(A, x, y)
    fails = fails + FSUNMatSpace_Band(A, lenrw, leniw)

    !======= Cleanup ===========
    call FSUNMatDestroy_Band(A)
    call FSUNMatDestroy_Band(B)
    call FN_VDestroy_Serial(x)
    call FN_VDestroy_Serial(y)

  end function smoke_tests

  integer(c_int) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding

    use fnvector_serial_mod
    use fsunmatrix_band_mod

    use test_sunmatrix

    implicit none

    type(SUNMatrix), pointer   :: A, I
    type(N_Vector), pointer   :: x, y
    real(c_double), pointer   :: Adata(:), Idata(:), xdata(:), ydata(:)
    integer(kind=myindextype) :: ii, jj, smu, istart, iend, offset

    fails = 0

    A => FSUNBandMatrix(N, mu, ml, sunctx)
    I => FSUNBandMatrix(N, 0_myindextype, 0_myindextype, sunctx)
    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)

    ! Fill identity matrix
    Idata => FSUNBandMatrix_Data(I)
    do jj = 1, N
      Idata(jj) = ONE
    end do

    ! Fill A matrix
    smu = FSUNBandMatrix_StoredUpperBandwidth(A)
    Adata => FSUNBandMatrix_Data(A)
    do jj = 1, N
      offset = (jj - 1)*(smu + ml + 1) + smu + 1 ! offset to diagonal
      istart = merge(-mu, -(jj - 1), jj > mu) ! above diagonal
      iend = merge(N - jj, ml, jj > N - ml)  ! below diagonal
      do ii = istart, iend
        Adata(offset + ii) = (jj - 1) - ii
      end do
    end do

    xdata => FN_VGetArrayPointer(x)
    ydata => FN_VGetArrayPointer(y)

    ! Fill vectors
    do jj = 0, N - 1
      ! x vector
      xdata(jj + 1) = jj
      ! y vector
      ydata(jj + 1) = ZERO
      istart = max(0_myindextype, jj - ml)
      iend = min(N - 1, jj + mu)
      do ii = istart, iend
        ydata(jj + 1) = ydata(jj + 1) + (ii + ii - jj)*(ii)
      end do
    end do

    fails = fails + Test_FSUNMatGetID(A, SUNMATRIX_BAND, 0)
    fails = fails + Test_FSUNMatClone(A, 0)
    fails = fails + Test_FSUNMatCopy(A, 0)
    fails = fails + Test_FSUNMatZero(A, 0)
    fails = fails + Test_FSUNMatScaleAdd(A, I, 0)
    fails = fails + Test_FSUNMatScaleAddI(A, I, 0)
    fails = fails + Test_FSUNMatMatvec(A, x, y, 0)
    fails = fails + Test_FSUNMatSpace(A, 0)

    ! cleanup
    call FSUNMatDestroy(A)
    call FSUNMatDestroy(I)
    call FN_VDestroy(x)
    call FN_VDestroy(y)

  end function unit_tests

end module

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunmatrix_band

  !======== Declarations ========
  implicit none
  integer(c_int) :: fails = 0

  !============== Introduction =============
  print *, 'Band SUNMatrix Fortran 2003 interface test'

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

! exported functions used by test_sunmatrix
integer(c_int) function check_matrix(B, A, tol) result(fails)
  use, intrinsic :: iso_c_binding

  use fsunmatrix_band_mod
  use test_utilities

  implicit none

  type(SUNMatrix) :: A, B
  real(c_double)  :: tol
  real(c_double), pointer :: Adata(:), Bdata(:)
  integer(c_long) :: N, smu, mu, ml, ii, istart, iend, jj, offset

  fails = 0

  N = FSUNBandMatrix_Columns(A)
  smu = FSUNBandMatrix_StoredUpperBandwidth(A)
  mu = FSUNBandMatrix_UpperBandwidth(A)
  ml = FSUNBandMatrix_LowerBandwidth(A)

  if (FSUNMatGetID(A) /= FSUNMatGetID(B)) then
    fails = 1
    return
  end if

  if (FSUNBandMatrix_Columns(A) /= FSUNBandMatrix_Columns(B)) then
    fails = 1
    return
  end if

  if (FSUNBandMatrix_Rows(A) /= FSUNBandMatrix_Rows(B)) then
    fails = 1
    return
  end if

  if (FSUNBandMatrix_LowerBandwidth(A) /= FSUNBandMatrix_LowerBandwidth(B)) then
    fails = 1
    return
  end if

  if (FSUNBandMatrix_UpperBandwidth(A) /= FSUNBandMatrix_UpperBandwidth(B)) then
    fails = 1
    return
  end if

  if (FSUNBandMatrix_StoredUpperBandwidth(A) /= FSUNBandMatrix_StoredUpperBandwidth(B)) then
    fails = 1
    return
  end if

  Adata => FSUNBandMatrix_Data(A)
  Bdata => FSUNBandMatrix_Data(B)
  do jj = 1, N
    offset = (jj - 1)*(smu + ml + 1) + smu + 1  ! offset to diagonal
    istart = merge(-mu, -(jj - 1), jj > mu) ! above diagonal
    iend = merge(N - jj, ml, jj > N - ml)  ! below diagonal
    do ii = istart, iend
      fails = fails + FNEQTOL(Adata(offset + ii), Bdata(offset + ii), tol)
    end do
  end do

end function check_matrix

integer(c_int) function check_matrix_entry(A, c, tol) result(fails)
  use, intrinsic :: iso_c_binding

  use fsunmatrix_band_mod
  use test_utilities

  implicit none

  type(SUNMatrix) :: A
  real(c_double)  :: c, tol
  real(c_double), pointer :: Adata(:)
  integer(c_long) :: N, smu, mu, ml, ii, istart, iend, jj, offset

  fails = 0

  N = FSUNBandMatrix_Columns(A)
  smu = FSUNBandMatrix_StoredUpperBandwidth(A)
  mu = FSUNBandMatrix_UpperBandwidth(A)
  ml = FSUNBandMatrix_LowerBandwidth(A)
  Adata => FSUNBandMatrix_Data(A)

  do jj = 1, N
    offset = (jj - 1)*(smu + ml + 1) + smu + 1  ! offset to diagonal
    istart = merge(-mu, -(jj - 1), jj > mu) ! above diagonal
    iend = merge(N - jj, ml, jj > N - ml)  ! below diagonal
    do ii = istart, iend
      if (FNEQTOL(Adata(offset + ii), c, tol) /= 0) then
        fails = fails + 1
        write (*, '(A,E10.1,A,E14.7,A,I9,A,E14.7)') "tol = ", tol, &
          "   c = ", c, "   data[", offset + ii, "] = ", Adata(offset + ii)
      end if
    end do
  end do

end function check_matrix_entry

logical function is_square(A) result(res)
  use fsundials_core_mod
  use fsunmatrix_band_mod

  implicit none

  type(SUNMatrix) :: A

  if (FSUNBandMatrix_Rows(A) == FSUNBandMatrix_Columns(A)) then
    res = .true.
  else
    res = .false.
  end if

end function is_square
