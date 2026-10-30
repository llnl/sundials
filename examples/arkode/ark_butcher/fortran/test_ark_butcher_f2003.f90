! ------------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
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
! Tests the ARKButcherTable F2003 interface.
! ------------------------------------------------------------------

module test_arkode_butcher_table

#if defined(SUNDIALS_INT32_T)
  integer, parameter :: myindextype = selected_int_kind(8)
#elif defined(SUNDIALS_INT64_T)
  integer, parameter :: myindextype = selected_int_kind(16)
#endif

contains

  integer function smoke_tests() result(ret)

    !======== Inclusions ==========
    use, intrinsic :: iso_c_binding
    use farkode_mod

    !======== Declarations ========
    implicit none
    type(c_ptr) :: ERK, DIRK
    integer(C_INT) :: ierr, q(1), p(1)
    integer(kind=myindextype) :: liw(1), lrw(1)
    real(C_DOUBLE) :: b(2), c(2), d(2), A(4)

    !===== Setup ====

    ! ARKODE_HEUN_EULER_2_1_2
    A = 0.d0
    b = 0.d0
    c = 0.d0
    d = 0.d0

    A(3) = 1.0d0
    b(1) = 0.5d0
    b(2) = 0.5d0
    c(2) = 1.0d0
    d(1) = 1.0d0

    !===== Test =====
    ERK = FARkodeButcherTable_LoadERK(ARKODE_HEUN_EULER_2_1_2)
    DIRK = FARkodeButcherTable_LoadDIRK(ARKODE_SDIRK_2_1_2)
    ierr = FARkodeButcherTable_CheckOrder(ERK, q, p, C_NULL_PTR)
    ierr = FARkodeButcherTable_CheckARKOrder(ERK, DIRK, q, p, C_NULL_PTR)
    call FARKodeButcherTable_Space(ERK, liw, lrw)
    call FARKodeButcherTable_Free(ERK)
    call FARKodeButcherTable_Free(DIRK)

    ERK = FARkodeButcherTable_Create(2, 2, 1, c, A, b, d)
    DIRK = FARkodeButcherTable_Alloc(2, 1)
    call FARKodeButcherTable_Free(DIRK)
    DIRK = FARkodeButcherTable_Copy(ERK)

    !==== Cleanup =====
    call FARKodeButcherTable_Free(ERK)
    call FARKodeButcherTable_Free(DIRK)

    ret = 0

  end function smoke_tests

end module

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_arkode_butcher_table

  !======== Declarations ========
  implicit none
  integer :: fails = 0

  !============== Introduction =============
  print *, 'ARKodeButcherTable Fortran 2003 interface test'

  fails = smoke_tests()
  if (fails /= 0) then
    print *, 'FAILURE: smoke tests failed'
    stop 1
  else
    print *, 'SUCCESS: smoke tests passed'
  end if

end program main
