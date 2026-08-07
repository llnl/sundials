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
! A helper module for the fortran tests
! -----------------------------------------------------------------

module test_utilities

  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
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

  real(c_double), parameter :: SUN_UNIT_ROUNDOFF = epsilon(1.0d0)

  real(c_double) :: NEG_TWO = -2.0d0
  real(c_double) :: NEG_ONE = -1.0d0
  real(c_double) :: NEG_HALF = -0.50d0
  real(c_double) :: ZERO = 0.0d0
  real(c_double) :: HALF = 0.5d0
  real(c_double) :: ONE = 1.0d0
  real(c_double) :: TWO = 2.0d0
  real(c_double) :: THREE = 3.0d0
  real(c_double) :: FOUR = 4.0d0
  real(c_double) :: FIVE = 5.0d0
  real(c_double) :: SIX = 6.0d0

  type(c_ptr)    :: sunctx

contains

  subroutine Test_Init(comm)
    implicit none
    integer(c_int), value :: comm
    integer(c_int)        :: retval

    retval = FSUNContext_Create(comm, sunctx)
    if (retval /= 0) then
      print *, 'ERROR in Test_Init: FSUNContext_Create returned nonzero'
      stop 1
    end if

  end subroutine

  subroutine Test_Finalize()
    implicit none
    integer(c_int) :: retval

    retval = FSUNContext_Free(sunctx)

  end subroutine

  integer(c_int) function FNEQTOL(a, b, tol) result(nequal)
    implicit none
    real(c_double) :: a, b, tol

    if (a /= a) then
      nequal = 1
    else if ((abs(a - b)/abs(b)) > tol) then
      nequal = 1
    else
      nequal = 0
    end if

  end function FNEQTOL

  integer(c_int) function FNEQ(a, b) result(nequal)
    implicit none
    real(c_double) :: a, b

    if (a /= a) then
      nequal = 1
    else if ((abs(a - b)/abs(b)) > (10*SUN_UNIT_ROUNDOFF)) then
      nequal = 1
    else
      nequal = 0
    end if
  end function FNEQ

end module
