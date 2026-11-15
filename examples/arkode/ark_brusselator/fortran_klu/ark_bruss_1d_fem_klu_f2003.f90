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
! Example problem:
!
! The following test simulates a brusselator problem from chemical
! kinetics.  This is a PDE system with 3 components, Y = [u,v,w],
! satisfying the equations,
!    du/dt = du*u_xx + a - (w+1)*u + v*u^2
!    dv/dt = dv*v_xx + w*u - v*u^2
!    dw/dt = dw*w_xx + (b-w)/ep - w*u
! for t in the interval [0, 10], x in [0, 10], with initial
! conditions
!    u(0,x) =  a  + 0.1*sin(pi*x),
!    v(0,x) = b/a + 0.1*sin(pi*x),
!    w(0,x) =  b  + 0.1*sin(pi*x),
! and with stationary boundary conditions, i.e.
!    u_t(t,0) = u_t(t,1) = 0,
!    v_t(t,0) = v_t(t,1) = 0,
!    w_t(t,0) = w_t(t,1) = 0.
!
! Here, we use a piecewise linear Galerkin finite element
! discretization in space, where all element-wise integrals are
! computed using 3-node Gaussian quadrature (since we will have
! quartic polynomials in the reaction terms for the u_t and v_t
! equations (including the test function)).  The time derivative
! terms in this system will include a mass matrix, giving rise to
! an ODE system of the form
!      M y_t = L y + R(y),
! where M is the 3x3 block mass matrix for each component, L is
! the 3x3 block Laplace operator for each component, and R(y) is
! comprised of the nonlinear reaction terms for each component.
! Since it it highly inefficient to rewrite this system as
!      y_t = M^{-1}(L y + R(y)),
! we solve this system using FARKODE, with a user-supplied mass
! matrix.  We therefore provide functions to evaluate the ODE RHS
!    f(t,y) = L y + R(y),
! its Jacobian
!    J(t,y) = L + dR/dy,
! and the mass matrix, M.
!
! We use N=201 spatial nodes, with parameters
!    a=0.6,  b=2.0,  du=0.025,  dv=0.025,  dw=0.025,  ep=1.d-5
!
! This program solves the problem with the DIRK method, using a
! Newton iteration with the SUNKLU sparse linear solvers for both
! the system and mass matrices.  These matrices are stored in
! compressed-sparse-row format.
!
! Output is printed 10 times throughout the defined time interval.
! Run statistics (optional outputs) are printed at the end.
! ------------------------------------------------------------------

module Bruss1DFEMKLU_UserData

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  ! number of equations
  integer(c_int64_t), parameter :: neqreal = 3

  ! ODE parameters
  integer(c_int64_t), parameter :: N = 201          ! number of intervals
  integer(c_int64_t), parameter :: neq = neqreal*N  ! set overall problem size
  integer(c_int64_t), parameter :: nnz = 15*neq
  real(c_double), parameter :: a = 0.6d0        ! constant forcing on u
  real(c_double), parameter :: b = 2.d0         ! steady-state value of w
  real(c_double), parameter :: du = 2.5d-2      ! diffusion coeff for u
  real(c_double), parameter :: dv = 2.5d-2      ! diffusion coeff for v
  real(c_double), parameter :: dw = 2.5d-2      ! diffusion coeff for w
  real(c_double), parameter :: ep = 1.d-5       ! stiffness parameter
  real(c_double), dimension(N) :: x             ! mesh node locations

contains

  ! function that maps 2D data into 1D address space
  ! (0-based since CSR matrix will be sent to C solver)
  integer(c_int64_t) function idx(ix, ivar)
    integer(c_int64_t):: ix
    integer(c_int) :: ivar
    idx = neqreal*(ix - 1) + ivar - 1
  end function idx

end module Bruss1DFEMKLU_UserData
! ------------------------------------------------------------------

! finite element basis functions
module FEMBasis
  use, intrinsic :: iso_c_binding

contains

  ! left/right basis functions
  real(c_double) function ChiL(xl, xr, x)
    real(c_double) :: xl, xr, x
    ChiL = (xr - x)/(xr - xl)
  end function ChiL

  real(c_double) function ChiR(xl, xr, x)
    real(c_double) :: xl, xr, x
    ChiR = (x - xl)/(xr - xl)
  end function ChiR

  ! derivatives of left/right basis functions
  real(c_double) function ChiL_x(xl, xr)
    real(c_double) :: xl, xr
    ChiL_x = 1.d0/(xl - xr)
  end function ChiL_X

  real(c_double) function ChiR_x(xl, xr)
    real(c_double) :: xl, xr
    ChiR_x = 1.d0/(xr - xl)
  end function ChiR_x

  ! FEM output evaluation routines: value and derivative
  real(c_double) function Eval(ul, ur, xl, xr, x)
    real(c_double) :: ul, ur, xl, xr, x
    Eval = ul*ChiL(xl, xr, x) + ur*ChiR(xl, xr, x)
  end function Eval

  real(c_double) function Eval_x(ul, ur, xl, xr)
    real(c_double) :: ul, ur, xl, xr
    Eval_x = ul*ChiL_x(xl, xr) + ur*ChiR_x(xl, xr)
  end function Eval_x

end module FEMBasis
! ------------------------------------------------------------------

! quadrature data
module Quadrature
  use, intrinsic :: iso_c_binding

contains

  ! nodes
  real(c_double) function X1(xl, xr)
    real(c_double) :: xl, xr
    X1 = 0.5d0*(xl + xr) - 0.5d0*(xr - xl)*0.774596669241483377035853079956d0
  end function X1

  real(c_double) function X2(xl, xr)
    real(c_double) :: xl, xr
    X2 = 0.5d0*(xl + xr)
  end function X2

  real(c_double) function X3(xl, xr)
    real(c_double) :: xl, xr
    X3 = 0.5d0*(xl + xr) + 0.5d0*(xr - xl)*0.774596669241483377035853079956d0
  end function X3

  ! quadrature
  real(c_double) function Quad(f1, f2, f3, xl, xr)
    real(c_double) :: f1, f2, f3, xl, xr
    real(c_double), parameter :: wt1 = 0.55555555555555555555555555555556d0
    real(c_double), parameter :: wt2 = 0.88888888888888888888888888888889d0
    real(c_double), parameter :: wt3 = 0.55555555555555555555555555555556d0
    Quad = 0.5d0*(xr - xl)*(wt1*f1 + wt2*f2 + wt3*f3)
  end function Quad

end module Quadrature
! ------------------------------------------------------------------

module bruss1D_ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use Bruss1DFEMKLU_UserData
  use fsundials_core_mod

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
  integer(c_int) function ImpRhsFn(tn, sunvec_y, sunvec_f, user_data) &
    result(ierr) bind(C)

    !======= Inclusions ===========
    use FEMBasis
    use Quadrature

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! Local data
    integer(c_int64_t) :: ix
    logical        :: left, right
    real(c_double) :: ul, ur, vl, vr, wl, wr, xl, xr, u, v, w, f1, f2, f3

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neqreal, N) :: yvec(:, :)
    real(c_double), pointer, dimension(neqreal, N) :: fvec(:, :)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_f)

    ! clear out rhs
    fvec = 0.d0

    ! iterate over intervals, filling in rhs function
    do ix = 1, N - 1

      ! set booleans to determine whether equations exist on the left/right */
      left = .true.
      right = .true.
      if (ix == 1) left = .false.
      if (ix == (N - 1)) right = .false.

      ! set nodal value shortcuts (interval index aligns with left node)
      ul = yvec(1, ix)
      vl = yvec(2, ix)
      wl = yvec(3, ix)
      ur = yvec(1, ix + 1)
      vr = yvec(2, ix + 1)
      wr = yvec(3, ix + 1)

      ! set mesh shortcuts
      xl = x(ix)
      xr = x(ix + 1)

      !    left test function
      if (left) then

        ! u -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl, xr))
        v = Eval(vl, vr, xl, xr, X1(xl, xr))
        w = Eval(wl, wr, xl, xr, X1(xl, xr))
        f1 = (a - (w + 1.d0)*u + v*u*u)*ChiL(xl, xr, X1(xl, xr))
        u = Eval(ul, ur, xl, xr, X2(xl, xr))
        v = Eval(vl, vr, xl, xr, X2(xl, xr))
        w = Eval(wl, wr, xl, xr, X2(xl, xr))
        f2 = (a - (w + 1.d0)*u + v*u*u)*ChiL(xl, xr, X2(xl, xr))
        u = Eval(ul, ur, xl, xr, X3(xl, xr))
        v = Eval(vl, vr, xl, xr, X3(xl, xr))
        w = Eval(wl, wr, xl, xr, X3(xl, xr))
        f3 = (a - (w + 1.d0)*u + v*u*u)*ChiL(xl, xr, X3(xl, xr))
        fvec(1, ix) = fvec(1, ix) + Quad(f1, f2, f3, xl, xr)

        ! u -- diffusion
        f1 = -du*Eval_x(ul, ur, xl, xr)*ChiL_x(xl, xr)
        fvec(1, ix) = fvec(1, ix) + Quad(f1, f1, f1, xl, xr)

        ! v -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl, xr))
        v = Eval(vl, vr, xl, xr, X1(xl, xr))
        w = Eval(wl, wr, xl, xr, X1(xl, xr))
        f1 = (w*u - v*u*u)*ChiL(xl, xr, X1(xl, xr))
        u = Eval(ul, ur, xl, xr, X2(xl, xr))
        v = Eval(vl, vr, xl, xr, X2(xl, xr))
        w = Eval(wl, wr, xl, xr, X2(xl, xr))
        f2 = (w*u - v*u*u)*ChiL(xl, xr, X2(xl, xr))
        u = Eval(ul, ur, xl, xr, X3(xl, xr))
        v = Eval(vl, vr, xl, xr, X3(xl, xr))
        w = Eval(wl, wr, xl, xr, X3(xl, xr))
        f3 = (w*u - v*u*u)*ChiL(xl, xr, X3(xl, xr))
        fvec(2, ix) = fvec(2, ix) + Quad(f1, f2, f3, xl, xr)

        ! v -- diffusion
        f1 = -dv*Eval_x(vl, vr, xl, xr)*ChiL_x(xl, xr)
        fvec(2, ix) = fvec(2, ix) + Quad(f1, f1, f1, xl, xr)

        ! w -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl, xr))
        v = Eval(vl, vr, xl, xr, X1(xl, xr))
        w = Eval(wl, wr, xl, xr, X1(xl, xr))
        f1 = ((b - w)/ep - w*u)*ChiL(xl, xr, X1(xl, xr))
        u = Eval(ul, ur, xl, xr, X2(xl, xr))
        v = Eval(vl, vr, xl, xr, X2(xl, xr))
        w = Eval(wl, wr, xl, xr, X2(xl, xr))
        f2 = ((b - w)/ep - w*u)*ChiL(xl, xr, X2(xl, xr))
        u = Eval(ul, ur, xl, xr, X3(xl, xr))
        v = Eval(vl, vr, xl, xr, X3(xl, xr))
        w = Eval(wl, wr, xl, xr, X3(xl, xr))
        f3 = ((b - w)/ep - w*u)*ChiL(xl, xr, X3(xl, xr))
        fvec(3, ix) = fvec(3, ix) + Quad(f1, f2, f3, xl, xr)

        ! w -- diffusion
        f1 = -dw*Eval_x(wl, wr, xl, xr)*ChiL_x(xl, xr)
        fvec(3, ix) = fvec(3, ix) + Quad(f1, f1, f1, xl, xr)

      end if

      !    right test function
      if (right) then

        ! u -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl, xr))
        v = Eval(vl, vr, xl, xr, X1(xl, xr))
        w = Eval(wl, wr, xl, xr, X1(xl, xr))
        f1 = (a - (w + 1.d0)*u + v*u*u)*ChiR(xl, xr, X1(xl, xr))
        u = Eval(ul, ur, xl, xr, X2(xl, xr))
        v = Eval(vl, vr, xl, xr, X2(xl, xr))
        w = Eval(wl, wr, xl, xr, X2(xl, xr))
        f2 = (a - (w + 1.d0)*u + v*u*u)*ChiR(xl, xr, X2(xl, xr))
        u = Eval(ul, ur, xl, xr, X3(xl, xr))
        v = Eval(vl, vr, xl, xr, X3(xl, xr))
        w = Eval(wl, wr, xl, xr, X3(xl, xr))
        f3 = (a - (w + 1.d0)*u + v*u*u)*ChiR(xl, xr, X3(xl, xr))
        fvec(1, ix + 1) = fvec(1, ix + 1) + Quad(f1, f2, f3, xl, xr)

        ! u -- diffusion
        f1 = -du*Eval_x(ul, ur, xl, xr)*ChiR_x(xl, xr)
        fvec(1, ix + 1) = fvec(1, ix + 1) + Quad(f1, f1, f1, xl, xr)

        ! v -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl, xr))
        v = Eval(vl, vr, xl, xr, X1(xl, xr))
        w = Eval(wl, wr, xl, xr, X1(xl, xr))
        f1 = (w*u - v*u*u)*ChiR(xl, xr, X1(xl, xr))
        u = Eval(ul, ur, xl, xr, X2(xl, xr))
        v = Eval(vl, vr, xl, xr, X2(xl, xr))
        w = Eval(wl, wr, xl, xr, X2(xl, xr))
        f2 = (w*u - v*u*u)*ChiR(xl, xr, X2(xl, xr))
        u = Eval(ul, ur, xl, xr, X3(xl, xr))
        v = Eval(vl, vr, xl, xr, X3(xl, xr))
        w = Eval(wl, wr, xl, xr, X3(xl, xr))
        f3 = (w*u - v*u*u)*ChiR(xl, xr, X3(xl, xr))
        fvec(2, ix + 1) = fvec(2, ix + 1) + Quad(f1, f2, f3, xl, xr)

        ! v -- diffusion
        f1 = -dv*Eval_x(vl, vr, xl, xr)*ChiR_x(xl, xr)
        fvec(2, ix + 1) = fvec(2, ix + 1) + Quad(f1, f1, f1, xl, xr)

        ! w -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl, xr))
        v = Eval(vl, vr, xl, xr, X1(xl, xr))
        w = Eval(wl, wr, xl, xr, X1(xl, xr))
        f1 = ((b - w)/ep - w*u)*ChiR(xl, xr, X1(xl, xr))
        u = Eval(ul, ur, xl, xr, X2(xl, xr))
        v = Eval(vl, vr, xl, xr, X2(xl, xr))
        w = Eval(wl, wr, xl, xr, X2(xl, xr))
        f2 = ((b - w)/ep - w*u)*ChiR(xl, xr, X2(xl, xr))
        u = Eval(ul, ur, xl, xr, X3(xl, xr))
        v = Eval(vl, vr, xl, xr, X3(xl, xr))
        w = Eval(wl, wr, xl, xr, X3(xl, xr))
        f3 = ((b - w)/ep - w*u)*ChiR(xl, xr, X3(xl, xr))
        fvec(3, ix + 1) = fvec(3, ix + 1) + Quad(f1, f2, f3, xl, xr)

        ! w -- diffusion
        f1 = -dw*Eval_x(wl, wr, xl, xr)*ChiR_x(xl, xr)
        fvec(3, ix + 1) = fvec(3, ix + 1) + Quad(f1, f1, f1, xl, xr)

      end if

    end do

    ! return success
    ierr = 0
    return

  end function ImpRhsFn
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Jac: The Jacobian function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function Jac(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                              sunvec_t1, sunvec_t2, sunvec_t3) result(ierr) bind(C, name='Jac')

    !======= Inclusions ===========
    use FEMBasis
    use Quadrature
    use fnvector_serial_mod
    use fsunmatrix_sparse_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(SUNMatrix)       :: sunmat_J  ! Jacobian SUNMatrix
    type(c_ptr), value :: user_data ! user-defined data
    type(N_Vector)        :: sunvec_t1 ! temporary N_Vectors
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! Local data
    integer(c_int64_t) :: ix, nz, Nint
    real(c_double) :: ul, uc, ur, vl, vc, vr, wl, wc, wr, xl, xc, xr
    real(c_double) :: u1, u2, u3, v1, v2, v3, w1, w2, w3
    real(c_double) :: df1, df2, df3, dQdf1, dQdf2, dQdf3
    real(c_double) :: ChiL1, ChiL2, ChiL3, ChiR1, ChiR2, ChiR3
    real(c_double), dimension(3, -1:1) :: Ju, Jv, Jw

    ! pointers to data in SUNDIALS vectors
    integer(c_int64_t), pointer, dimension(nnz)       :: Jcolvals(:)
    integer(c_int64_t), pointer, dimension(neq + 1)     :: Jrowptrs(:)
    real(c_double), pointer, dimension(nnz)       :: Jdata(:)
    real(c_double), pointer, dimension(neqreal, N) :: yvec(:, :)
    real(c_double), pointer, dimension(neqreal, N) :: fvec(:, :)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_f)
    Jdata(1:nnz) => FSUNSparseMatrix_Data(sunmat_J)
    Jcolvals(1:nnz) => FSUNSparseMatrix_IndexValues(sunmat_J)
    Jrowptrs(1:neq + 1) => FSUNSparseMatrix_IndexPointers(sunmat_J)

    ! check that vector/matrix dimensions match up
    if ((3*N /= neq) .or. (nnz < 27*(N - 2))) then
      ierr = 1
      return
    end if

    ! set integer*4 version of N for call to idx()
    Nint = N

    ! clear out Jacobian matrix data
    Jdata = 0.d0
    nz = 0

    ! Dirichlet boundary at left
    Jrowptrs(idx(1_c_int64_t, 1) + 1) = nz
    Jrowptrs(idx(1_c_int64_t, 2) + 1) = nz
    Jrowptrs(idx(1_c_int64_t, 3) + 1) = nz

    ! iterate through nodes, filling in matrix by rows
    do ix = 2, N - 1

      ! set nodal value shortcuts (interval index aligns with left node)
      xl = x(ix - 1)
      ul = yvec(1, ix - 1)
      vl = yvec(2, ix - 1)
      wl = yvec(3, ix - 1)
      xc = x(ix)
      uc = yvec(1, ix)
      vc = yvec(2, ix)
      wc = yvec(3, ix)
      xr = x(ix + 1)
      ur = yvec(1, ix + 1)
      vr = yvec(2, ix + 1)
      wr = yvec(3, ix + 1)

      ! compute entries of all Jacobian rows at node ix
      Ju = 0.d0
      Jv = 0.d0
      Jw = 0.d0

      ! first compute dependence on values to left and center

      !    evaluate relevant variables in left subinterval
      u1 = Eval(ul, uc, xl, xc, X1(xl, xc))
      v1 = Eval(vl, vc, xl, xc, X1(xl, xc))
      w1 = Eval(wl, wc, xl, xc, X1(xl, xc))
      u2 = Eval(ul, uc, xl, xc, X2(xl, xc))
      v2 = Eval(vl, vc, xl, xc, X2(xl, xc))
      w2 = Eval(wl, wc, xl, xc, X2(xl, xc))
      u3 = Eval(ul, uc, xl, xc, X3(xl, xc))
      v3 = Eval(vl, vc, xl, xc, X3(xl, xc))
      w3 = Eval(wl, wc, xl, xc, X3(xl, xc))

      dQdf1 = Quad(1.d0, 0.d0, 0.d0, xl, xc)
      dQdf2 = Quad(0.d0, 1.d0, 0.d0, xl, xc)
      dQdf3 = Quad(0.d0, 0.d0, 1.d0, xl, xc)

      ChiL1 = ChiL(xl, xc, X1(xl, xc))
      ChiL2 = ChiL(xl, xc, X2(xl, xc))
      ChiL3 = ChiL(xl, xc, X3(xl, xc))
      ChiR1 = ChiR(xl, xc, X1(xl, xc))
      ChiR2 = ChiR(xl, xc, X2(xl, xc))
      ChiR3 = ChiR(xl, xc, X3(xl, xc))

      !    compute diffusion Jacobian components

      !    L_u = -du * u_x * ChiR_x
      !     dL_u/dul
      Ju(1, -1) = (-du)*Quad(1.d0, 1.d0, 1.d0, xl, xc)*ChiL_x(xl, xc)*ChiR_x(xl, xc)
      !     dL_u/duc
      Ju(1, 0) = (-du)*Quad(1.d0, 1.d0, 1.d0, xl, xc)*ChiR_x(xl, xc)*ChiR_x(xl, xc)

      !    L_v = -dv * v_x * ChiR_x
      !     dL_v/dvl
      Jv(2, -1) = (-dv)*Quad(1.d0, 1.d0, 1.d0, xl, xc)*ChiL_x(xl, xc)*ChiR_x(xl, xc)
      !     dL_v/dvc
      Jv(2, 0) = (-dv)*Quad(1.d0, 1.d0, 1.d0, xl, xc)*ChiR_x(xl, xc)*ChiR_x(xl, xc)

      !    L_w =  -dw * w_x * ChiR_x
      !     dL_w/dwl
      Jw(3, -1) = (-dw)*Quad(1.d0, 1.d0, 1.d0, xl, xc)*ChiL_x(xl, xc)*ChiR_x(xl, xc)
      !     dL_w/dwc
      Jw(3, 0) = (-dw)*Quad(1.d0, 1.d0, 1.d0, xl, xc)*ChiR_x(xl, xc)*ChiR_x(xl, xc)

      !    compute reaction Jacobian components

      !    R_u = (a - (w+1.d0)*u + v*u*u)
      !     dR_u/dul
      df1 = (-(w1 + 1.d0) + 2.d0*v1*u1)*ChiL1*ChiR1
      df2 = (-(w2 + 1.d0) + 2.d0*v2*u2)*ChiL2*ChiR2
      df3 = (-(w3 + 1.d0) + 2.d0*v3*u3)*ChiL3*ChiR3
      Ju(1, -1) = Ju(1, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/duc
      df1 = (-(w1 + 1.d0) + 2.d0*v1*u1)*ChiR1*ChiR1
      df2 = (-(w2 + 1.d0) + 2.d0*v2*u2)*ChiR2*ChiR2
      df3 = (-(w3 + 1.d0) + 2.d0*v3*u3)*ChiR3*ChiR3
      Ju(1, 0) = Ju(1, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dvl
      df1 = (u1*u1)*ChiL1*ChiR1
      df2 = (u2*u2)*ChiL2*ChiR2
      df3 = (u3*u3)*ChiL3*ChiR3
      Ju(2, -1) = Ju(2, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dvc
      df1 = (u1*u1)*ChiR1*ChiR1
      df2 = (u2*u2)*ChiR2*ChiR2
      df3 = (u3*u3)*ChiR3*ChiR3
      Ju(2, 0) = Ju(2, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dwl
      df1 = (-u1)*ChiL1*ChiR1
      df2 = (-u2)*ChiL2*ChiR2
      df3 = (-u3)*ChiL3*ChiR3
      Ju(3, -1) = Ju(3, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dwc
      df1 = (-u1)*ChiR1*ChiR1
      df2 = (-u2)*ChiR2*ChiR2
      df3 = (-u3)*ChiR3*ChiR3
      Ju(3, 0) = Ju(3, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !    R_v = (w*u - v*u*u)
      !     dR_v/dul
      df1 = (w1 - 2.d0*v1*u1)*ChiL1*ChiR1
      df2 = (w2 - 2.d0*v2*u2)*ChiL2*ChiR2
      df3 = (w3 - 2.d0*v3*u3)*ChiL3*ChiR3
      Jv(1, -1) = Jv(1, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/duc
      df1 = (w1 - 2.d0*v1*u1)*ChiR1*ChiR1
      df2 = (w2 - 2.d0*v2*u2)*ChiR2*ChiR2
      df3 = (w3 - 2.d0*v3*u3)*ChiR3*ChiR3
      Jv(1, 0) = Jv(1, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dvl
      df1 = (-u1*u1)*ChiL1*ChiR1
      df2 = (-u2*u2)*ChiL2*ChiR2
      df3 = (-u3*u3)*ChiL3*ChiR3
      Jv(2, -1) = Jv(2, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dvc
      df1 = (-u1*u1)*ChiR1*ChiR1
      df2 = (-u2*u2)*ChiR2*ChiR2
      df3 = (-u3*u3)*ChiR3*ChiR3
      Jv(2, 0) = Jv(2, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dwl
      df1 = (u1)*ChiL1*ChiR1
      df2 = (u2)*ChiL2*ChiR2
      df3 = (u3)*ChiL3*ChiR3
      Jv(3, -1) = Jv(3, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dwc
      df1 = (u1)*ChiR1*ChiR1
      df2 = (u2)*ChiR2*ChiR2
      df3 = (u3)*ChiR3*ChiR3
      Jv(3, 0) = Jv(3, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !    R_w = ((b-w)/ep - w*u)
      !     dR_w/dul
      df1 = (-w1)*ChiL1*ChiR1
      df2 = (-w2)*ChiL2*ChiR2
      df3 = (-w3)*ChiL3*ChiR3
      Jw(1, -1) = Jw(1, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_w/duc
      df1 = (-w1)*ChiR1*ChiR1
      df2 = (-w2)*ChiR2*ChiR2
      df3 = (-w3)*ChiR3*ChiR3
      Jw(1, 0) = Jw(1, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_w/dwl
      df1 = (-1.d0/ep - u1)*ChiL1*ChiR1
      df2 = (-1.d0/ep - u2)*ChiL2*ChiR2
      df3 = (-1.d0/ep - u3)*ChiL3*ChiR3
      Jw(3, -1) = Jw(3, -1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_w/dwc
      df1 = (-1.d0/ep - u1)*ChiR1*ChiR1
      df2 = (-1.d0/ep - u2)*ChiR2*ChiR2
      df3 = (-1.d0/ep - u3)*ChiR3*ChiR3
      Jw(3, 0) = Jw(3, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      ! second compute dependence on values to center and right

      !    evaluate relevant variables in right subinterval
      u1 = Eval(uc, ur, xc, xr, X1(xc, xr))
      v1 = Eval(vc, vr, xc, xr, X1(xc, xr))
      w1 = Eval(wc, wr, xc, xr, X1(xc, xr))
      u2 = Eval(uc, ur, xc, xr, X2(xc, xr))
      v2 = Eval(vc, vr, xc, xr, X2(xc, xr))
      w2 = Eval(wc, wr, xc, xr, X2(xc, xr))
      u3 = Eval(uc, ur, xc, xr, X3(xc, xr))
      v3 = Eval(vc, vr, xc, xr, X3(xc, xr))
      w3 = Eval(wc, wr, xc, xr, X3(xc, xr))

      dQdf1 = Quad(1.d0, 0.d0, 0.d0, xc, xr)
      dQdf2 = Quad(0.d0, 1.d0, 0.d0, xc, xr)
      dQdf3 = Quad(0.d0, 0.d0, 1.d0, xc, xr)

      ChiL1 = ChiL(xc, xr, X1(xc, xr))
      ChiL2 = ChiL(xc, xr, X2(xc, xr))
      ChiL3 = ChiL(xc, xr, X3(xc, xr))
      ChiR1 = ChiR(xc, xr, X1(xc, xr))
      ChiR2 = ChiR(xc, xr, X2(xc, xr))
      ChiR3 = ChiR(xc, xr, X3(xc, xr))

      !    compute diffusion Jacobian components

      !    L_u = -du * u_x * ChiL_x
      !     dL_u/duc
      Ju(1, 0) = Ju(1, 0) + (-du)*Quad(1.d0, 1.d0, 1.d0, xc, xr)*ChiL_x(xc, xr)*ChiL_x(xc, xr)

      !     dL_u/dur
      Ju(1, 1) = Ju(1, 1) + (-du)*Quad(1.d0, 1.d0, 1.d0, xc, xr)*ChiL_x(xc, xr)*ChiR_x(xc, xr)

      !    L_v = -dv * v_x * ChiL_x
      !     dL_v/dvc
      Jv(2, 0) = Jv(2, 0) + (-dv)*Quad(1.d0, 1.d0, 1.d0, xc, xr)*ChiL_x(xc, xr)*ChiL_x(xc, xr)

      !     dL_v/dvr
      Jv(2, 1) = Jv(2, 1) + (-dv)*Quad(1.d0, 1.d0, 1.d0, xc, xr)*ChiL_x(xc, xr)*ChiR_x(xc, xr)

      !    L_w =  -dw * w_x * ChiL_x
      !     dL_w/dwc
      Jw(3, 0) = Jw(3, 0) + (-dw)*Quad(1.d0, 1.d0, 1.d0, xc, xr)*ChiL_x(xc, xr)*ChiL_x(xc, xr)

      !     dL_w/dwr
      Jw(3, 1) = Jw(3, 1) + (-dw)*Quad(1.d0, 1.d0, 1.d0, xc, xr)*ChiL_x(xc, xr)*ChiR_x(xc, xr)

      !    compute reaction Jacobian components

      !    R_u = (a - (w+1.d0)*u + v*u*u)
      !     dR_u/duc
      df1 = (-(w1 + 1.d0) + 2.d0*v1*u1)*ChiL1*ChiL1
      df2 = (-(w2 + 1.d0) + 2.d0*v2*u2)*ChiL2*ChiL2
      df3 = (-(w3 + 1.d0) + 2.d0*v3*u3)*ChiL3*ChiL3
      Ju(1, 0) = Ju(1, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dur
      df1 = (-(w1 + 1.d0) + 2.d0*v1*u1)*ChiL1*ChiR1
      df2 = (-(w2 + 1.d0) + 2.d0*v2*u2)*ChiL2*ChiR2
      df3 = (-(w3 + 1.d0) + 2.d0*v3*u3)*ChiL3*ChiR3
      Ju(1, 1) = Ju(1, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dvc
      df1 = (u1*u1)*ChiL1*ChiL1
      df2 = (u2*u2)*ChiL2*ChiL2
      df3 = (u3*u3)*ChiL3*ChiL3
      Ju(2, 0) = Ju(2, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dvr
      df1 = (u1*u1)*ChiL1*ChiR1
      df2 = (u2*u2)*ChiL2*ChiR2
      df3 = (u3*u3)*ChiL3*ChiR3
      Ju(2, 1) = Ju(2, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dwc
      df1 = (-u1)*ChiL1*ChiL1
      df2 = (-u2)*ChiL2*ChiL2
      df3 = (-u3)*ChiL3*ChiL3
      Ju(3, 0) = Ju(3, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_u/dwr
      df1 = (-u1)*ChiL1*ChiR1
      df2 = (-u2)*ChiL2*ChiR2
      df3 = (-u3)*ChiL3*ChiR3
      Ju(3, 1) = Ju(3, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !    R_v = (w*u - v*u*u)
      !     dR_v/duc
      df1 = (w1 - 2.d0*v1*u1)*ChiL1*ChiL1
      df2 = (w2 - 2.d0*v2*u2)*ChiL2*ChiL2
      df3 = (w3 - 2.d0*v3*u3)*ChiL3*ChiL3
      Jv(1, 0) = Jv(1, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dur
      df1 = (w1 - 2.d0*v1*u1)*ChiL1*ChiR1
      df2 = (w2 - 2.d0*v2*u2)*ChiL2*ChiR2
      df3 = (w3 - 2.d0*v3*u3)*ChiL3*ChiR3
      Jv(1, 1) = Jv(1, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dvc
      df1 = (-u1*u1)*ChiL1*ChiL1
      df2 = (-u2*u2)*ChiL2*ChiL2
      df3 = (-u3*u3)*ChiL3*ChiL3
      Jv(2, 0) = Jv(2, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dvr
      df1 = (-u1*u1)*ChiL1*ChiR1
      df2 = (-u2*u2)*ChiL2*ChiR2
      df3 = (-u3*u3)*ChiL3*ChiR3
      Jv(2, 1) = Jv(2, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dwc
      df1 = (u1)*ChiL1*ChiL1
      df2 = (u2)*ChiL2*ChiL2
      df3 = (u3)*ChiL3*ChiL3
      Jv(3, 0) = Jv(3, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_v/dwr
      df1 = (u1)*ChiL1*ChiR1
      df2 = (u2)*ChiL2*ChiR2
      df3 = (u3)*ChiL3*ChiR3
      Jv(3, 1) = Jv(3, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !    R_w = ((b-w)/ep - w*u)
      !     dR_w/duc
      df1 = (-w1)*ChiL1*ChiL1
      df2 = (-w2)*ChiL2*ChiL2
      df3 = (-w3)*ChiL3*ChiL3
      Jw(1, 0) = Jw(1, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_w/dur
      df1 = (-w1)*ChiL1*ChiR1
      df2 = (-w2)*ChiL2*ChiR2
      df3 = (-w3)*ChiL3*ChiR3
      Jw(1, 1) = Jw(1, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_w/dwc
      df1 = (-1.d0/ep - u1)*ChiL1*ChiL1
      df2 = (-1.d0/ep - u2)*ChiL2*ChiL2
      df3 = (-1.d0/ep - u3)*ChiL3*ChiL3
      Jw(3, 0) = Jw(3, 0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      !     dR_w/dwr
      df1 = (-1.d0/ep - u1)*ChiL1*ChiR1
      df2 = (-1.d0/ep - u2)*ChiL2*ChiR2
      df3 = (-1.d0/ep - u3)*ChiL3*ChiR3
      Jw(3, 1) = Jw(3, 1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

      ! insert Jacobian entries into CSR matrix structure

      !   Ju row
      Jrowptrs(idx(ix, 1) + 1) = nz

      Jdata(nz + 1:nz + 3) = (/Ju(1, -1), Ju(2, -1), Ju(3, -1)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix - 1, 1), idx(ix - 1, 2), idx(ix - 1, 3)/)
      nz = nz + 3

      Jdata(nz + 1:nz + 3) = (/Ju(1, 0), Ju(2, 0), Ju(3, 0)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix, 1), idx(ix, 2), idx(ix, 3)/)
      nz = nz + 3

      Jdata(nz + 1:nz + 3) = (/Ju(1, 1), Ju(2, 1), Ju(3, 1)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix + 1, 1), idx(ix + 1, 2), idx(ix + 1, 3)/)
      nz = nz + 3

      !   Jv row
      Jrowptrs(idx(ix, 2) + 1) = nz

      Jdata(nz + 1:nz + 3) = (/Jv(1, -1), Jv(2, -1), Jv(3, -1)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix - 1, 1), idx(ix - 1, 2), idx(ix - 1, 3)/)
      nz = nz + 3

      Jdata(nz + 1:nz + 3) = (/Jv(1, 0), Jv(2, 0), Jv(3, 0)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix, 1), idx(ix, 2), idx(ix, 3)/)
      nz = nz + 3

      Jdata(nz + 1:nz + 3) = (/Jv(1, 1), Jv(2, 1), Jv(3, 1)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix + 1, 1), idx(ix + 1, 2), idx(ix + 1, 3)/)
      nz = nz + 3

      !   Jw row
      Jrowptrs(idx(ix, 3) + 1) = nz

      Jdata(nz + 1:nz + 3) = (/Jw(1, -1), Jw(2, -1), Jw(3, -1)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix - 1, 1), idx(ix - 1, 2), idx(ix - 1, 3)/)
      nz = nz + 3

      Jdata(nz + 1:nz + 3) = (/Jw(1, 0), Jw(2, 0), Jw(3, 0)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix, 1), idx(ix, 2), idx(ix, 3)/)
      nz = nz + 3

      Jdata(nz + 1:nz + 3) = (/Jw(1, 1), Jw(2, 1), Jw(3, 1)/)
      Jcolvals(nz + 1:nz + 3) = (/idx(ix + 1, 1), idx(ix + 1, 2), idx(ix + 1, 3)/)
      nz = nz + 3

    end do

    ! Dirichlet boundary at right
    Jrowptrs(idx(Nint, 1) + 1) = nz
    Jrowptrs(idx(Nint, 2) + 1) = nz
    Jrowptrs(idx(Nint, 3) + 1) = nz

    ! signal end of data in CSR matrix
    Jrowptrs(idx(Nint, 3) + 2) = nz

    ! return success
    ierr = 0
    return

  end function Jac
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Mass matrix computation routine
  ! ----------------------------------------------------------------
  integer(c_int) function Mass(tn, sunmat_M, user_data, &
                               sunvec_t1, sunvec_t2, sunvec_t3) result(ierr) bind(C, name='Mass')

    !======= Inclusions ===========
    use FEMBasis
    use Quadrature
    use fnvector_serial_mod
    use fsunmatrix_sparse_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn         ! current time
    type(SUNMatrix)       :: sunmat_M   ! Jacobian SUNMatrix
    type(c_ptr), value :: user_data  ! user-defined data
    type(N_Vector)        :: sunvec_t1  ! temporary N_Vectors
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! Local data
    integer(c_int64_t) :: ix, nz, Nint
    real(c_double) :: xl, xc, xr, Ml, Mc, Mr, ChiL1, ChiL2, ChiL3, ChiR1, ChiR2, ChiR3
    logical        :: left, right

    ! pointers to data in SUNDIALS vectors
    integer(c_int64_t), pointer, dimension(nnz)   :: Mcolvals(:)
    integer(c_int64_t), pointer, dimension(neq + 1) :: Mrowptrs(:)
    real(c_double), pointer, dimension(nnz)   :: Mdata(:)

    !======= Internals ============
    ! get data arrays from SUNDIALS vectors
    Mdata(1:nnz) => FSUNSparseMatrix_Data(sunmat_M)
    Mcolvals(1:nnz) => FSUNSparseMatrix_IndexValues(sunmat_M)
    Mrowptrs(1:neq + 1) => FSUNSparseMatrix_IndexPointers(sunmat_M)

    ! check that vector/matrix dimensions match up
    if ((3*N /= neq) .or. (nnz /= 15*neq)) then
      ierr = 1
      return
    end if

    ! set integer*4 version of N for call to idx()
    Nint = N

    ! clear out Jacobian matrix data
    Mdata = 0.d0
    nz = 0

    ! iterate through nodes, filling in matrix by rows
    do ix = 1, N

      ! set booleans to determine whether intervals exist on the left/right */
      left = .true.
      right = .true.
      if (ix == 1) left = .false.
      if (ix == N) right = .false.

      ! set nodal value shortcuts (interval index aligns with left node)
      if (left) then
        xl = x(ix - 1)
      end if
      xc = x(ix)
      if (right) then
        xr = x(ix + 1)
      end if

      ! compute entries of all mass matrix rows at node ix
      Ml = 0.d0
      Mc = 0.d0
      Mr = 0.d0

      ! first compute dependence on values to left and center
      if (left) then

        ChiL1 = ChiL(xl, xc, X1(xl, xc))
        ChiL2 = ChiL(xl, xc, X2(xl, xc))
        ChiL3 = ChiL(xl, xc, X3(xl, xc))
        ChiR1 = ChiR(xl, xc, X1(xl, xc))
        ChiR2 = ChiR(xl, xc, X2(xl, xc))
        ChiR3 = ChiR(xl, xc, X3(xl, xc))

        Ml = Ml + Quad(ChiL1*ChiR1, ChiL2*ChiR2, ChiL3*ChiR3, xl, xc)
        Mc = Mc + Quad(ChiR1*ChiR1, ChiR2*ChiR2, ChiR3*ChiR3, xl, xc)

      end if

      ! second compute dependence on values to center and right
      if (right) then

        ChiL1 = ChiL(xc, xr, X1(xc, xr))
        ChiL2 = ChiL(xc, xr, X2(xc, xr))
        ChiL3 = ChiL(xc, xr, X3(xc, xr))
        ChiR1 = ChiR(xc, xr, X1(xc, xr))
        ChiR2 = ChiR(xc, xr, X2(xc, xr))
        ChiR3 = ChiR(xc, xr, X3(xc, xr))

        Mc = Mc + Quad(ChiL1*ChiL1, ChiL2*ChiL2, ChiL3*ChiL3, xc, xr)
        Mr = Mr + Quad(ChiL1*ChiR1, ChiL2*ChiR2, ChiL3*ChiR3, xc, xr)

      end if

      ! insert mass matrix entries into CSR matrix structure

      !   u row
      Mrowptrs(idx(ix, 1) + 1) = nz
      if (left) then
        nz = nz + 1
        Mdata(nz) = Ml
        Mcolvals(nz) = idx(ix - 1, 1)
      end if
      nz = nz + 1
      Mdata(nz) = Mc
      Mcolvals(nz) = idx(ix, 1)
      if (right) then
        nz = nz + 1
        Mdata(nz) = Mr
        Mcolvals(nz) = idx(ix + 1, 1)
      end if

      !   v row
      Mrowptrs(idx(ix, 2) + 1) = nz
      if (left) then
        nz = nz + 1
        Mdata(nz) = Ml
        Mcolvals(nz) = idx(ix - 1, 2)
      end if
      nz = nz + 1
      Mdata(nz) = Mc
      Mcolvals(nz) = idx(ix, 2)
      if (right) then
        nz = nz + 1
        Mdata(nz) = Mr
        Mcolvals(nz) = idx(ix + 1, 2)
      end if

      !   w row
      Mrowptrs(idx(ix, 3) + 1) = nz
      if (left) then
        nz = nz + 1
        Mdata(nz) = Ml
        Mcolvals(nz) = idx(ix - 1, 3)
      end if
      nz = nz + 1
      Mdata(nz) = Mc
      Mcolvals(nz) = idx(ix, 3)
      if (right) then
        nz = nz + 1
        Mdata(nz) = Mr
        Mcolvals(nz) = idx(ix + 1, 3)
      end if

    end do

    ! signal end of data in CSR matrix
    Mrowptrs(idx(Nint, 3) + 2) = nz

    ! return success
    ierr = 0
    return

  end function Mass
  ! ----------------------------------------------------------------

end module bruss1D_ode_mod
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use bruss1D_ode_mod            ! custom problem-specification module
  use farkode_mod                ! Fortran interface to the ARKODE module
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep module
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_sparse_mod      ! Fortran interface to sparse SUNMatrix
  use fsunlinsol_klu_mod         ! Fortran interface to dense SUNLinearSolver
  use Bruss1DFEMKLU_UserData     ! Declarations and indexing

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)     :: ctx         ! SUNDIALS context for the simulation
  real(c_double)  :: tstart      ! initial time
  real(c_double)  :: tend        ! final time
  real(c_double)  :: rtol, atol  ! relative and absolute tolerance
  real(c_double)  :: dtout       ! output time interval
  real(c_double)  :: tout        ! output time
  real(c_double)  :: tcur(1)     ! current time
  real(c_double)  :: pi, h, z    ! constants and variables to help with mesh
  integer(c_int)  :: time_dep    ! time dependence parameter
  integer(c_int)  :: ierr        ! error flag from C functions
  integer(c_int)  :: nout        ! number of outputs
  integer(c_int)  :: outstep     ! output loop counter
  integer(c_int)  :: sparsetype  ! CSR signal, here
  integer(c_long) :: mxsteps     ! max num steps
  integer(c_int64_t) :: i

  type(N_Vector), pointer :: sunvec_y    ! sundials vector
  type(N_Vector), pointer :: sunvec_u    ! sundials vector
  type(N_Vector), pointer :: sunvec_v    ! sundials vector
  type(N_Vector), pointer :: sunvec_w    ! sundials vector
  type(SUNMatrix), pointer :: sunmat_A    ! sundials (linsol) matrix
  type(SUNMatrix), pointer :: sunmat_M    ! sundials (mass) matrix
  type(SUNLinearSolver), pointer :: sunls_A     ! sundials linear solver
  type(SUNLinearSolver), pointer :: sunls_M     ! sundials linear solver
  type(c_ptr)                    :: arkode_mem  ! ARKODE memory
  type(c_ptr)                    :: outstr      ! standard output file stream
  real(c_double), pointer, dimension(neqreal, N) :: yvec(:, :)   ! underlying vector y
  real(c_double), pointer, dimension(neqreal, N) :: umask(:, :)  ! identifier for u
  real(c_double), pointer, dimension(neqreal, N) :: vmask(:, :)  ! identifier for v
  real(c_double), pointer, dimension(neqreal, N) :: wmask(:, :)  ! identifier for w

  !======= Internals ============

  ! create the SUNDIALS context
  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)

  ! initialize ODE
  tstart = 0.0d0
  tend = 10.0d0
  tcur = tstart
  tout = tstart
  dtout = (tend - tstart)/10.d0
  nout = ceiling(tend/dtout)

  ! create and assign SUNDIALS N_Vectors
  sunvec_y => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_y)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  yvec(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_y)

  sunvec_u => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_u)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  umask(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_u)

  sunvec_v => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_v)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  vmask(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_v)

  sunvec_w => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_w)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  wmask(1:neqreal, 1:N) => FN_VGetArrayPointer(sunvec_w)

  ! set up spatial mesh; this [arbitrarily] clusters
  ! more intervals near the end points of the interval
  pi = 4.d0*atan(1.d0)
  h = 10.d0/(N - 1)
  do i = 1, N
    z = -5.d0 + h*(i - 1)
    x(i) = 0.5d0/atan(5.d0)*atan(z) + 0.5d0
  end do

  ! output mesh to disk
  open (200, file='bruss_FEM_mesh.txt')
  do i = 1, N
    write (200, *) x(i)
  end do
  close (200)

  ! set initial conditions into yvec
  do i = 1, N
    yvec(1, i) = a + 0.1d0*sin(pi*x(i))   ! u0
    yvec(2, i) = b/a + 0.1d0*sin(pi*x(i))   ! v0
    yvec(3, i) = b + 0.1d0*sin(pi*x(i))   ! w0
  end do

  ! set mask values for each solution component
  umask = 0.d0
  vmask = 0.d0
  wmask = 0.d0
  do i = 1, N
    umask(1, i) = 1.d0
    vmask(2, i) = 1.d0
    wmask(3, i) = 1.d0
  end do

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(ImpRhsFn), tstart, sunvec_y, ctx)
  if (.not. c_associated(arkode_mem)) print *, 'ERROR: arkode_mem = NULL'

  ! Tell ARKODE to use a sparse linear solver for both Newton and mass matrix systems.
  sparsetype = 1
  sunmat_A => FSUNSparseMatrix(neq, neq, nnz, sparsetype, ctx)
  if (.not. associated(sunmat_A)) then
    print *, 'ERROR: sunmat_A = NULL'
    stop 1
  end if

  sunmat_M => FSUNSparseMatrix(neq, neq, nnz, sparsetype, ctx)
  if (.not. associated(sunmat_M)) then
    print *, 'ERROR: sunmat_M = NULL'
    stop 1
  end if

  sunls_A => FSUNLinSol_KLU(sunvec_y, sunmat_A, ctx)
  if (.not. associated(sunls_A)) then
    print *, 'ERROR: sunls_A = NULL'
    stop 1
  end if

  ierr = FARKodeSetLinearSolver(arkode_mem, sunls_A, sunmat_A)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetLinearSolver'
    stop 1
  end if

  ierr = FARKodeSetJacFn(arkode_mem, c_funloc(Jac))
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetJacFn'
    stop 1
  end if

  sunls_M => FSUNLinSol_KLU(sunvec_y, sunmat_M, ctx)
  if (.not. associated(sunls_M)) then
    print *, 'ERROR: sunls_M = NULL'
    stop 1
  end if

  time_dep = 0
  ierr = FARKodeSetMassLinearSolver(arkode_mem, sunls_M, sunmat_M, time_dep)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetMassLinearSolver'
    stop 1
  end if

  ierr = FARKodeSetMassFn(arkode_mem, c_funloc(Mass))
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetMassFn'
    stop 1
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-6
  atol = 1.0d-11

  ierr = FARKodeSStolerances(arkode_mem, rtol, atol)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSStolerances, ierr = ', ierr, '; halting'
    stop 1
  end if

  ! set residual tolerance with the same atol as above
  ierr = FARKodeResStolerance(arkode_mem, atol)
  if (ierr /= 0) then
    print *, 'Error in FARKodeResStolerance, ierr = ', ierr, '; halting'
    stop 1
  end if

  ! Set maximum number of internal time steps
  mxsteps = 1000
  ierr = FARKodeSetMaxNumSteps(arkode_mem, mxsteps)
  if (ierr /= 0) then
    print *, 'Error in FARKodeSetNonlinConvCoef'
    stop 1
  end if

  ! Open output stream for results
  open (501, file='bruss_FEM_u.txt')
  open (502, file='bruss_FEM_v.txt')
  open (503, file='bruss_FEM_w.txt')

  ! output initial condition to disk
  write (501, *) (yvec(1, i), i=1, N)
  write (502, *) (yvec(2, i), i=1, N)
  write (503, *) (yvec(3, i), i=1, N)

  ! output solver parameters to screen
  ierr = FSUNDIALSFileOpen('stdout', 'w', outstr)
  if (ierr /= 0) then
    print *, 'Error in FSUNDIALSFileOpen'
    stop 1
  end if
  ierr = FARKodeWriteParameters(arkode_mem, outstr)
  if (ierr /= 0) then
    print *, 'Error in FARKodeWriteParameters'
    stop 1
  end if
  ierr = FSUNDIALSFileClose(outstr)
  if (ierr /= 0) then
    print *, 'Error in FSUNDIALSFileClose'
    stop 1
  end if

  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '        t         ||u||_rms    ||v||_rms    ||w||_rms'
  print *, '  ----------------------------------------------------'
  print '(3x,4(es12.5,1x))', tcur, sqrt(sum(yvec*yvec*umask)/N), &
    sqrt(sum(yvec*yvec*vmask)/N), sqrt(sum(yvec*yvec*wmask)/N)
  do outstep = 1, nout

    ! set the next output time
    tout = min(tout + dtout, tend)

    ierr = FARKodeSetStopTime(arkode_mem, tout)
    if (ierr /= 0) then
      print *, 'Error in FARKodeSetStopTime, ierr = ', ierr, '; halting'
      stop 1
    end if

    ! call ARKodeEvolve
    ierr = FARKodeEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
    if (ierr < 0) then
      print *, 'Error in FARKodeEvolve, ierr = ', ierr, '; halting'
      stop 1
    end if

    ! output current solution information (using yvec)
    print '(3x,4(es12.5,1x))', Tcur, sqrt(sum(yvec*yvec*umask)/N), &
      sqrt(sum(yvec*yvec*vmask)/N), sqrt(sum(yvec*yvec*wmask)/N)
    write (501, *) (yvec(1, i), i=1, N)
    write (502, *) (yvec(2, i), i=1, N)
    write (503, *) (yvec(3, i), i=1, N)

  end do
  print *, '  ----------------------------------------------------'
  close (501)
  close (502)
  close (503)

  ! diagnostics output
  call ARKStepStats(arkode_mem)

  ! clean up
  call FARKodeFree(arkode_mem)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_v)
  call FN_VDestroy(sunvec_w)
  call FSUNMatDestroy(sunmat_A)
  call FSUNMatDestroy(sunmat_M)
  ierr = FSUNLinSolFree(sunls_A)
  ierr = FSUNLinSolFree(sunls_M)
  ierr = FSUNContext_Free(ctx)

end program main
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! ARKStepStats
!
! Print ARKODE statstics to stdandard out
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
  integer(c_long) :: nlinsetups(1) ! num linear solver setups
  integer(c_long) :: netfails(1)   ! num error test fails

  real(c_double)  :: hlast(1)      ! last step size
  real(c_double)  :: hcur(1)       ! step size for next step
  real(c_double)  :: tcur(1)       ! internal time reached

  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nncfails(1)   ! nonlinear solver fails
  integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations
  integer(c_long) :: nmassevals(1) ! number of Mass matrix evaluations

  !======= Internals ============

  ierr = FARKodeGetNumSteps(arkode_mem, nsteps)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumSteps, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumStepAttempts, retval = ', ierr, '; halting'
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

  ierr = FARKodeGetLastStep(arkode_mem, hlast)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetLastStep, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetCurrentStep(arkode_mem, hcur)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetCurrentStep, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetCurrentTime(arkode_mem, tcur)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetCurrentTime, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumLinSolvSetups, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumErrTestFails(arkode_mem, netfails)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumErrTestFails, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumNonlinSolvIters(arkode_mem, nniters)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvIters, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumNonlinSolvConvFails, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumJacEvals(arkode_mem, njacevals)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumJacEvals, retval = ', ierr, '; halting'
    stop 1
  end if

  ierr = FARKodeGetNumMassSetups(arkode_mem, nmassevals)
  if (ierr /= 0) then
    print *, 'Error in FARKodeGetNumMassSetups, retval = ', ierr, '; halting'
    stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)', 'Total internal steps taken    =', nsteps
  print '(4x,A,i9)', 'Total internal steps attempts =', nst_a
  print '(4x,A,i9)', 'Total rhs exp function calls  =', nfe
  print '(4x,A,i9)', 'Total rhs imp function calls  =', nfi
  print '(4x,A,i9)', 'Total jac function calls      =', njacevals
  print '(4x,A,i9)', 'Total mass function calls     =', nmassevals
  print '(4x,A,i9)', 'Num lin solver setup calls    =', nlinsetups
  print '(4x,A,i9)', 'Num error test failures       =', netfails
  print '(4x,A,es12.5)', 'Last internal step size       =', hlast
  print '(4x,A,es12.5)', 'Next internal step size       =', hcur
  print '(4x,A,es12.5)', 'Current internal time         =', tcur
  print '(4x,A,i9)', 'Num nonlinear solver iters    =', nniters
  print '(4x,A,i9)', 'Num nonlinear solver fails    =', nncfails
  print *, ' '

  return

end subroutine ARKStepStats
! ----------------------------------------------------------------
