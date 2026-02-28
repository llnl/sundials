! ------------------------------------------------------------------
! Programmer(s): Steven Roberts @ LLNL
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
! Routine to test loading Butcher tables via strings
! ------------------------------------------------------------------

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  use farkode_mod

  !======= Declarations =========
  implicit none

  type(c_ptr)        :: table      ! Butcher table object
  character(len=19) :: table_name ! table name
  integer(c_int)     :: ierr       ! error flag
  integer(c_int)     :: q(1)       ! table order
  integer(c_int)     :: p(1)       ! table embedded order

  print *, 'Loading ARKODE_DORMAND_PRINCE_7_4_5'
  table = FARKodeButcherTable_LoadERKByName('ARKODE_DORMAND_PRINCE_7_4_5')
  if (.not. c_associated(table)) then
    write (error_unit, *) 'FARKodeButcherTable_LoadERKByName returned NULL'
    stop 1
  end if

  print *, 'Checking ARKODE_DORMAND_PRINCE_7_4_5 order'
  ierr = FARKodeButcherTable_CheckOrder(table, q, p, c_null_ptr); 
  if (ierr /= 0) then
    write (error_unit, *) 'FARKodeButcherTable_CheckOrder returned ', ierr
    stop 1
  end if

  call FARKodeButcherTable_Free(table)

  print *, 'Loading ARKODE_TRBDF2_3_3_2'
  table_name = 'ARKODE_TRBDF2_3_3_2'
  table = FARKodeButcherTable_LoadDIRKByName(table_name)
  if (.not. c_associated(table)) then
    write (error_unit, *) 'FARKodeButcherTable_LoadDIRKByName returned NULL'
    stop 1
  end if

  print *, 'Checking ARKODE_TRBDF2_3_3_2 order'
  ierr = FARKodeButcherTable_CheckOrder(table, q, p, c_null_ptr); 
  if (ierr /= 0) then
    write (error_unit, *) 'FARKodeButcherTable_CheckOrder returned ', ierr
    stop 1
  end if

  call FARKodeButcherTable_Free(table)

  print *, 'Loading ARKODE_DIRK_NONE'
  table = FARKodeButcherTable_LoadDIRKByName('ARKODE_DIRK_NONE')
  if (c_associated(table)) then
    write (error_unit, *) 'FARKodeButcherTable_LoadDIRKByName returned non-NULL for ARKODE_DIRK_NONE'
    stop 1
  end if

  print *, 'Loading invalid table. This should print an error'
  table = FARKodeButcherTable_LoadERKByName('does not exist')
  if (c_associated(table)) then
    write (error_unit, *) 'FARKodeButcherTable_LoadERKByName returned non-NULL for invalid table name'
    stop 1
  end if

  print *, 'SUCCESS'
end program main

