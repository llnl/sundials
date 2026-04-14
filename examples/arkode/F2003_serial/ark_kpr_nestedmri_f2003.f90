! ------------------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!
! Based on the ark_kpr_nestedmri.cpp in examples/arkode/CXX_serial
! ------------------------------------------------------------------------------
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
! ------------------------------------------------------------------------------
! Nested multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
!
!    [u]' = [ G   e   e ] [ (u^2 - p(t) - 2) / (2u) ] +  [ p'(t)/(2u) ]
!    [v]    [ e  al  be ] [ (v^2 - q(t) - 2) / (2v) ]    [ q'(t)/(2v) ]
!    [w]    [ e -be  al ] [ (w^2 - r(t) - 2) / (2w) ]    [ r'(t)/(2w) ]
!
! where
!
!    p(t) = 0.5 * cos(t),
!    q(t) = cos(om * t * (1 + exp(-(t - 2)^2))), and
!    r(t) = cos(om * om * t * (1 + exp(-(t - 3)^2))).
!
! The first row corresponds to the slowest time scale, the second row to the
! intermediate time scale, and third row the fast time scale.
!
! This problem has the analytical solution:
!
!   u(t) = sqrt(2 + p(t)),
!   v(t) = sqrt(2 + q(t)), and
!   w(t) = sqrt(2 + r(t)).
!
! The problem parameters are
!
!   G  = -10.0 : stiffness at slow time scale
!   e  =   0.5 : fast/slow coupling strength
!   al =  -1.0 : oscillatory coupling between v and w
!   be =   1.0 : oscillatory coupling between v and w
!   om =  50.0 : time-scale separation factor
!
! The stiffness of the slow time scale is determined by G and for |G| > 50 it is
! 'stiff' and ideally suited to a multirate method that is implicit at the slow
! time scale.
!
! The coupling strength between the slow and faster components is proportional
! to |e|. The coupling between the intermediate and fast components is
! determined by al and be.
!
! The "intermediate" variable, v, oscillates at a frequency "om" times faster
! than u, and the "fast" variable, w, oscillates at a frequency "om" times
! faster than v.
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! Module containing ODE parameters and right hand side functions
! ------------------------------------------------------------------------------

module kpr_mod

  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  implicit none

  ! Problem parameters
  real(c_double), parameter :: G = -10.0d0
  real(c_double), parameter :: e = 0.5d0
  real(c_double), parameter :: al = -1.0d0
  real(c_double), parameter :: be = 1.0d0
  real(c_double), parameter :: om = 50.0d0

contains

  ! ff routine to compute the fast portion of the ODE RHS
  integer(c_int) function ff(t, y, ydot, user_data) result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value :: t
    type(N_Vector)        :: y
    type(N_Vector)        :: ydot
    type(c_ptr), value    :: user_data

    real(c_double), pointer, dimension(3) :: ydata(:)
    real(c_double), pointer, dimension(3) :: ydotdata(:)
    real(c_double) :: u, v, w, tmp1, tmp2, tmp3

    ! Get vector data arrays
    ydata => FN_VGetArrayPointer(y)
    ydotdata => FN_VGetArrayPointer(ydot)

    ! Extract variables
    u = ydata(1)
    v = ydata(2)
    w = ydata(3)

    ! fill in the RHS function:
    !   [ 0   0   0 ] [ (u^2 - p - 2) / (2u) ] +  [      0         ]
    !   [ 0   0   0 ] [ (v^2 - q - 2) / (2v) ]    [      0         ]
    !   [ e -be  al ] [ (w^2 - r - 2) / (2w) ]    [ rdot(t) / (2w) ]
    tmp1 = (-2.0d0 + u*u - p(t))/(2.0d0*u)
    tmp2 = (-2.0d0 + v*v - q(t))/(2.0d0*v)
    tmp3 = (-2.0d0 + w*w - r(t))/(2.0d0*w)

    ydotdata(1) = 0.0d0
    ydotdata(2) = 0.0d0
    ydotdata(3) = e*tmp1 - be*tmp2 + al*tmp3 + rdot(t)/(2.0d0*w)

    ! Return success
    ierr = 0
    return

  end function ff

  ! fm routine to compute the intermediate portion of the ODE RHS
  integer(c_int) function fm(t, y, ydot, user_data) result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value :: t
    type(N_Vector)        :: y
    type(N_Vector)        :: ydot
    type(c_ptr), value    :: user_data

    real(c_double), pointer, dimension(3) :: ydata(:)
    real(c_double), pointer, dimension(3) :: ydotdata(:)
    real(c_double) :: u, v, w, tmp1, tmp2, tmp3

    ! Get vector data arrays
    ydata => FN_VGetArrayPointer(y)
    ydotdata => FN_VGetArrayPointer(ydot)

    ! Extract variables
    u = ydata(1)
    v = ydata(2)
    w = ydata(3)

    ! fill in the RHS function:
    !   [ 0   0   0 ] [ (u^2 - p - 2) / (2u) ] +  [      0         ]
    !   [ e  al  be ] [ (v^2 - q - 2) / (2v) ]    [ qdot(t) / (2v) ]
    !   [ 0   0   0 ] [ (w^2 - r - 2) / (2w) ]    [      0         ]
    tmp1 = (-2.0d0 + u*u - p(t))/(2.0d0*u)
    tmp2 = (-2.0d0 + v*v - q(t))/(2.0d0*v)
    tmp3 = (-2.0d0 + w*w - r(t))/(2.0d0*w)

    ydotdata(1) = 0.0d0
    ydotdata(2) = e*tmp1 + al*tmp2 + be*tmp3 + qdot(t)/(2.0d0*v)
    ydotdata(3) = 0.0d0

    ierr = 0
    return

  end function fm

  ! fs routine to compute the slow portion of the ODE RHS
  integer(c_int) function fs(t, y, ydot, user_data) result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value :: t
    type(N_Vector)        :: y
    type(N_Vector)        :: ydot
    type(c_ptr), value    :: user_data

    real(c_double), pointer, dimension(3) :: ydata(:)
    real(c_double), pointer, dimension(3) :: ydotdata(:)
    real(c_double) :: u, v, w, tmp1, tmp2, tmp3

    ! Get vector data arrays
    ydata => FN_VGetArrayPointer(y)
    ydotdata => FN_VGetArrayPointer(ydot)

    ! Extract variables
    u = ydata(1)
    v = ydata(2)
    w = ydata(3)

    ! fill in the RHS function:
    !   [ G   e   e ] [ (u^2 - p - 2) / (2u) ] +  [ pdot(t) / (2u) ]
    !   [ 0   0   0 ] [ (v^2 - q - 2) / (2v) ]    [      0         ]
    !   [ 0   0   0 ] [ (w^2 - r - 2) / (2w) ]    [      0         ]
    tmp1 = (-2.0d0 + u*u - p(t))/(2.0d0*u)
    tmp2 = (-2.0d0 + v*v - q(t))/(2.0d0*v)
    tmp3 = (-2.0d0 + w*w - r(t))/(2.0d0*w)

    ydotdata(1) = G*tmp1 + e*tmp2 + e*tmp3 + pdot(t)/(2.0d0*u)
    ydotdata(2) = 0.0d0
    ydotdata(3) = 0.0d0

    ierr = 0
    return

  end function fs

  ! ----------------
  ! Helper functions
  ! ----------------

  real(c_double) function p(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = 0.5d0*cos(t)
    return
  end function p

  real(c_double) function q(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = cos(om*t*(1.0d0 + exp(-(t - 2.0d0)*(t - 2.0d0))))
    return
  end function q

  real(c_double) function r(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = cos(om*om*t*(1.0d0 + exp(-(t - 3.0d0)*(t - 3.0d0))))
    return
  end function r

  real(c_double) function pdot(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = -0.5d0*sin(t)
    return
  end function pdot

  real(c_double) function qdot(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    real(c_double) :: tm2
    real(c_double) :: eterm
    tm2 = t - 2.0d0
    eterm = exp(-tm2*tm2)
    result = -sin(om*t*(1.0d0 + eterm))*om*(1.0d0 + eterm*(1.0d0 - 2.0d0*t*tm2))
    return
  end function qdot

  real(c_double) function rdot(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    real(c_double) :: tm3
    real(c_double) :: eterm
    tm3 = t - 3.0d0
    eterm = exp(-tm3*tm3)
    result = -sin(om*om*t*(1.0d0 + eterm))*om*om*(1.0d0 + eterm*(1.0d0 - 2.0d0*t*tm3))
  end function rdot

  real(c_double) function utrue(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = sqrt(2.0d0 + p(t))
    return
  end function utrue

  real(c_double) function vtrue(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = sqrt(2.0d0 + q(t))
    return
  end function vtrue

  real(c_double) function wtrue(t) result(result)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    result = sqrt(2.0d0 + r(t))
    return
  end function wtrue

  integer(c_int) function Ytrue(t, y) result(ierr)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t
    type(N_Vector) :: y
    real(c_double), pointer, dimension(3) :: ydata(:)
    ydata => FN_VGetArrayPointer(y)
    ydata(1) = utrue(t)
    ydata(2) = vtrue(t)
    ydata(3) = wtrue(t)
    ierr = 0
    return
  end function Ytrue

end module kpr_mod

! ------------------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------------------

program main

  use, intrinsic :: iso_fortran_env
  use farkode_mod
  use farkode_arkstep_mod
  use farkode_erkstep_mod
  use farkode_mristep_mod
  use fnvector_serial_mod
  use fsunadaptcontroller_soderlind_mod
  use fsunadaptcontroller_mrihtol_mod
  use kpr_mod
  implicit none

  type(c_ptr) sunctx ! SUNDIALS simulation context

  type(N_Vector), pointer :: yvec                   ! solution vector
  real(c_double), pointer, dimension(3) :: ydata(:) ! solution vector array

  ! ARKODE integrators
  type(c_ptr) :: f_arkode_mem = c_null_ptr ! fast ARKODE integrator
  type(c_ptr) :: m_arkode_mem = c_null_ptr ! intermediate ARKODE integrator
  type(c_ptr) :: s_arkode_mem = c_null_ptr ! slow ARKODE integrator

  ! ARKODE inner steppers
  type(c_ptr) :: f_stepper = c_null_ptr ! fast stepper
  type(c_ptr) :: m_stepper = c_null_ptr ! intermediate stepper

  ! Slow and intermediate error controllers
  type(SUNAdaptController), pointer :: m_controller_H   ! step controller
  type(SUNAdaptController), pointer :: m_controller_Tol ! tolerance controller
  type(SUNAdaptController), pointer :: m_controller     ! h-tol controller

  type(SUNAdaptController), pointer :: s_controller_H   ! step controller
  type(SUNAdaptController), pointer :: s_controller_Tol ! tolerance controller
  type(SUNAdaptController), pointer :: s_controller     ! h-tol controller

  ! General problem parameters
#if defined(SUNDIALS_INT32_T)
  integer(kind=selected_int_kind(8)), parameter :: neq = 3
#elif defined(SUNDIALS_INT64_T)
  integer(kind=selected_int_kind(16)), parameter :: neq = 3
#endif
  real(c_double), parameter :: t0 = 0.0d0              ! initial time
  real(c_double), parameter :: tf = 5.0d0              ! final time
  integer(c_int), parameter :: num_output = 20                  ! number of outputs
  real(c_double), parameter :: rtol = 1.0d-4             ! relative tolerance
  real(c_double), parameter :: atol = 1.0d-11            ! absolute tolerance
  integer(c_int), parameter :: acc_type = ARK_ACCUMERROR_MAX ! error accumulation type

  ! Output variables
  real(c_double) :: output_dt        ! output interval
  real(c_double) :: tout             ! output time
  real(c_double) :: tret(1)          ! return time
  integer(c_int) :: iout             ! output counter
  real(c_double) :: uerr, verr, werr ! solution error

  ! Error return flag
  integer(c_int) :: retval = 0

  ! File pointer for output
  type(c_ptr) :: fp

  ! -----------------------
  ! Create SUNDIALS context
  ! -----------------------

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
  call check_retval(retval, 'FSUNContext_Create')

  ! -------------------------------
  ! Create initial condition vector
  ! -------------------------------

  yvec => FN_VNew_Serial(neq, sunctx)
  if (.not. associated(yvec)) then
    print *, 'ERROR: N_VNew_Serial failed'
    stop 1
  end if

  retval = Ytrue(T0, yvec)
  call check_retval(retval, "Ytrue")

  ! ----------------------
  ! Create fast integrator
  ! ----------------------

  f_arkode_mem = FERKStepCreate(c_funloc(ff), t0, yvec, sunctx)
  if (.not. c_associated(f_arkode_mem)) then
    print *, 'ERROR: f_arkode_mem = NULL'
    stop 1
  end if

  retval = FARKodeSStolerances(f_arkode_mem, rtol, atol)
  call check_retval(retval, "FARKodeSStolerances")

  retval = FARKodeSetOrder(f_arkode_mem, 4)
  call check_retval(retval, "FARKodeSetOrder")

  retval = FARKodeSetAccumulatedErrorType(f_arkode_mem, acc_type)
  call check_retval(retval, "FARKodeSetAccumulatedErrorType")

  retval = FARKodeCreateMRIStepInnerStepper(f_arkode_mem, f_stepper)
  call check_retval(retval, "FARKodeCreateMRIStepInnerStepper")

  ! ------------------------------
  ! Create intermediate integrator
  ! ------------------------------

  m_arkode_mem = FMRIStepCreate(c_funloc(fm), c_null_funptr, t0, yvec, f_stepper, sunctx)
  if (.not. c_associated(m_arkode_mem)) then
    print *, 'ERROR: m_arkode_mem = NULL'
    stop 1
  end if

  retval = FARKodeSStolerances(m_arkode_mem, rtol, atol)
  call check_retval(retval, "FARKodeSStolerances")

  retval = FARKodeSetOrder(m_arkode_mem, 4)
  call check_retval(retval, "FARKodeSetOrder")

  retval = FARKodeSetAccumulatedErrorType(m_arkode_mem, acc_type)
  call check_retval(retval, "FARKodeSetAccumulatedErrorType")

  m_controller_H => FSUNAdaptController_I(sunctx)
  if (.not. associated(m_controller_H)) then
    print *, 'ERROR: m_controller_H = NULL'
    stop 1
  end if

  m_controller_Tol => FSUNAdaptController_I(sunctx)
  if (.not. associated(m_controller_Tol)) then
    print *, 'ERROR: m_controller_Tol = NULL'
    stop 1
  end if

  m_controller => FSUNAdaptController_MRIHTol(m_controller_H, m_controller_Tol, sunctx)
  if (.not. associated(m_controller)) then
    print *, 'ERROR: m_controller = NULL'
    stop 1
  end if

  retval = FARKodeSetAdaptController(m_arkode_mem, m_controller)
  call check_retval(retval, "FARKodeSetAdaptController")

  retval = FARKodeCreateMRIStepInnerStepper(m_arkode_mem, m_stepper)
  call check_retval(retval, "FARKodeCreateMRIStepInnerStepper")

  ! ----------------------
  ! Create slow integrator
  ! ----------------------

  s_arkode_mem = FMRIStepCreate(c_funloc(fs), c_null_funptr, t0, yvec, m_stepper, sunctx)
  if (.not. c_associated(s_arkode_mem)) then
    print *, 'ERROR: s_arkode_mem = NULL'
    stop 1
  end if

  retval = FARKodeSStolerances(s_arkode_mem, rtol, atol)
  call check_retval(retval, "FARKodeSStolerances")

  retval = FARKodeSetOrder(s_arkode_mem, 4)
  call check_retval(retval, "FARKodeSetOrder")

  s_controller_H => FSUNAdaptController_I(sunctx)
  if (.not. associated(s_controller_H)) then
    print *, 'ERROR: s_controller_H = NULL'
    stop 1
  end if

  s_controller_Tol => FSUNAdaptController_I(sunctx)
  if (.not. associated(s_controller_Tol)) then
    print *, 'ERROR: s_controller_Tol = NULL'
    stop 1
  end if

  s_controller => FSUNAdaptController_MRIHTol(s_controller_H, s_controller_Tol, sunctx)
  if (.not. associated(s_controller)) then
    print *, 'ERROR: s_controller = NULL'
    stop 1
  end if

  retval = FARKodeSetAdaptController(s_arkode_mem, s_controller)
  call check_retval(retval, "FARKodeSetAdaptController")

  ! ---------------
  ! Advance in time
  ! ---------------

  output_dt = (tf - t0)/num_output ! output interval
  tout = t0 + output_dt         ! output time
  tret = t0                     ! return time

  ! Get state data to output
  ydata => FN_VGetArrayPointer(yvec)

  print '(7A23)', "t", "u", "v", "w", "u err", "v err", "w err"
  do iout = 1, 7*23
    write (*, '(A)', advance='no') "-"
  end do
  write (*, *)
  print '(7ES23.15)', tret(1), ydata(1), ydata(2), ydata(3), 0.0d0, 0.0d0, 0.0d0

  do iout = 1, num_output

    ! Stop at output time (do not interpolate)
    retval = FARKodeSetStopTime(s_arkode_mem, tout)

    retval = FARKodeEvolve(s_arkode_mem, tout, yvec, tret, ARK_NORMAL)
    call check_retval(retval, "FARKodeEvolve")

    ! Output solution and error
    uerr = abs(utrue(tret(1)) - ydata(1))
    verr = abs(vtrue(tret(1)) - ydata(2))
    werr = abs(wtrue(tret(1)) - ydata(3))

    print '(7ES23.15)', tret(1), ydata(1), ydata(2), ydata(3), uerr, verr, werr

    ! Update output time
    tout = tout + output_dt

  end do

  do iout = 1, 7*23
    write (*, '(A)', advance='no') "-"
  end do
  write (*, *)
  write (*, *)

  ! ----------------
  ! Final statistics
  ! ----------------

  ! Open up the file output.log for writing
  retval = FSUNDIALSFileOpen("stdout", "w+", fp)
  call check_retval(retval, "FSUNDIALSFileOpen")

  print '(A)', "Slow Integrator Stats"
  flush (output_unit)
  retval = FARKodePrintAllStats(s_arkode_mem, fp, SUN_OUTPUTFORMAT_TABLE)
  write (*, *)
  print '(A)', "Intermediate Integrator Stats"
  flush (output_unit)
  retval = FARKodePrintAllStats(m_arkode_mem, fp, SUN_OUTPUTFORMAT_TABLE)
  write (*, *)
  print '(A)', "Fast Integrator Stats"
  flush (output_unit)
  retval = FARKodePrintAllStats(f_arkode_mem, fp, SUN_OUTPUTFORMAT_TABLE)

  ! Close the file
  retval = FSUNDIALSFileClose(fp)
  call check_retval(retval, "FSUNDIALSFileOpen")

  ! --------
  ! Clean up
  ! --------

  ! Solution vector
  call FN_VDestroy(yvec)

  ! Fast integrator
  retval = FMRIStepInnerStepper_Free(f_stepper)
  call FARKodeFree(f_arkode_mem)

  ! Intermediate integrator
  retval = FSUNAdaptController_Destroy(m_controller)
  retval = FSUNAdaptController_Destroy(m_controller_H)
  retval = FSUNAdaptController_Destroy(m_controller_Tol)
  retval = FMRIStepInnerStepper_Free(m_stepper)
  call FARKodeFree(m_arkode_mem)

  ! Slow integrator
  retval = FSUNAdaptController_Destroy(s_controller)
  retval = FSUNAdaptController_Destroy(s_controller_H)
  retval = FSUNAdaptController_Destroy(s_controller_Tol)
  call FARKodeFree(s_arkode_mem)

  ! Context
  retval = FSUNContext_Free(sunctx)

end program main

! ------------------------------------------------------------------------------
! Utility to check return values
! ------------------------------------------------------------------------------

subroutine check_retval(retval, name)
  use, intrinsic :: iso_c_binding

  character(len=*) :: name
  integer(c_int)   :: retval

  if (retval < 0) then
    write (*, '(A,A,A,I3)') 'ERROR: ', name, ' returned ', retval
    stop 1
  end if
end subroutine check_retval
