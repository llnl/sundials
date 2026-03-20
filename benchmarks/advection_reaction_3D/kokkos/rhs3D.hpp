/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                David J. Gardner, Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------*/

#ifndef ADVECTION_REACTION_3D_RHS_HPP
#define ADVECTION_REACTION_3D_RHS_HPP

#include "advection_reaction_3D.hpp"

/* --------------------------------------------------------------
 * Right hand side (RHS) and residual functions
 * --------------------------------------------------------------*/

/* Compute the advection term f(t,y) = -c (grad * y). This is done using
   upwind 1st order finite differences.  At present, only periodic boundary
   conditions are supported, which are handled via MPI's Cartesian
   communicator (even for serial runs). */
static int Advection(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* set variable shortcuts */
  const int nxl        = udata->grid->nxl;
  const int nyl        = udata->grid->nyl;
  const int nzl        = udata->grid->nzl;
  const int dof        = udata->grid->dof;
  const sunrealtype c  = udata->c;
  const sunrealtype cx = -c / udata->grid->dx;
  const sunrealtype cy = -c / udata->grid->dy;
  const sunrealtype cz = -c / udata->grid->dz;

  /* local variables */
  int retval;

  /* fill send buffers and begin exchanging boundary information */
  SUNDIALS_MARK_BEGIN(udata->prof, "Neighbor Exchange");
  retval = FillSendBuffers(y, udata);
  if (check_retval(&retval, "FillSendBuffers", 1, udata->myid)) { return (-1); }
  retval = udata->grid->ExchangeStart();
  if (check_retval(&retval, "ExchangeStart", 1, udata->myid)) { return (-1); }
  SUNDIALS_MARK_END(udata->prof, "Neighbor Exchange");

  /* set output to zero */
  N_VConst(0.0, ydot);

  /* create 4D views of the state and RHS vectors */
  Vec4D Yview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y)), nxl, nyl,
              nzl, dof);
  Vec4D dYview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(ydot)), nxl,
               nyl, nzl, dof);

  /* iterate over domain interior, computing advection */
  if (c > 0.0)
  {
    /* flow moving in the positive x,y,z direction */
    Kokkos::parallel_for(
      "AdvectionInteriorRight", Range3D({1, 1, 1}, {nxl, nyl, nzl}),
      KOKKOS_LAMBDA(int i, int j, int k) {
        const sunrealtype u_ijk = Yview(i, j, k, 0);
        const sunrealtype v_ijk = Yview(i, j, k, 1);
        const sunrealtype w_ijk = Yview(i, j, k, 2);

        // grad * u
        dYview(i, j, k, 0) = cz * (u_ijk - Yview(i, j, k - 1, 0));  // du/dz
        dYview(i, j, k, 0) += cy * (u_ijk - Yview(i, j - 1, k, 0)); // du/dy
        dYview(i, j, k, 0) += cx * (u_ijk - Yview(i - 1, j, k, 0)); // du/dx

        // grad * v
        dYview(i, j, k, 1) = cz * (v_ijk - Yview(i, j, k - 1, 1));  // dv/dz
        dYview(i, j, k, 1) += cy * (v_ijk - Yview(i, j - 1, k, 1)); // dv/dy
        dYview(i, j, k, 1) += cx * (v_ijk - Yview(i - 1, j, k, 1)); // dv/dx

        // grad * w
        dYview(i, j, k, 2) = cz * (w_ijk - Yview(i, j, k - 1, 2));  // dw/dz
        dYview(i, j, k, 2) += cy * (w_ijk - Yview(i, j - 1, k, 2)); // dw/dy
        dYview(i, j, k, 2) += cx * (w_ijk - Yview(i - 1, j, k, 2)); // dw/dx
      });
  }
  else if (c < 0.0)
  {
    /* flow moving in the negative x,y,z direction */
    Kokkos::parallel_for(
      "AdvectionInteriorLeft", Range3D({0, 0, 0}, {nxl - 1, nyl - 1, nzl - 1}),
      KOKKOS_LAMBDA(int i, int j, int k) {
        const sunrealtype u_ijk = Yview(i, j, k, 0);
        const sunrealtype v_ijk = Yview(i, j, k, 1);
        const sunrealtype w_ijk = Yview(i, j, k, 2);

        // grad * u
        dYview(i, j, k, 0) = cz * (Yview(i, j, k + 1, 0) - u_ijk);  // du/dz
        dYview(i, j, k, 0) += cy * (Yview(i, j + 1, k, 0) - u_ijk); // du/dy
        dYview(i, j, k, 0) += cx * (Yview(i + 1, j, k, 0) - u_ijk); // du/dx

        // grad * v
        dYview(i, j, k, 1) = cz * (Yview(i, j, k + 1, 1) - v_ijk);  // dv/dz
        dYview(i, j, k, 1) += cy * (Yview(i, j + 1, k, 1) - v_ijk); // dv/dy
        dYview(i, j, k, 1) += cx * (Yview(i + 1, j, k, 1) - v_ijk); // dv/dx

        // grad * w
        dYview(i, j, k, 2) = cz * (Yview(i, j, k + 1, 2) - w_ijk);  // dw/dz
        dYview(i, j, k, 2) += cy * (Yview(i, j + 1, k, 2) - w_ijk); // dw/dy
        dYview(i, j, k, 2) += cx * (Yview(i + 1, j, k, 2) - w_ijk); // dw/dx
      });
  }

  /* finish exchanging boundary information */
  SUNDIALS_MARK_BEGIN(udata->prof, "Neighbor Exchange");
  retval = udata->grid->ExchangeEnd();
  if (check_retval(&retval, "ExchangeEnd", 1, udata->myid)) { return (-1); }
  SUNDIALS_MARK_END(udata->prof, "Neighbor Exchange");

  /* compute advection at process boundaries */
  if (c > 0.0)
  {
    /* Flow moving in the positive x,y,z direction:
       boundaries are west face, south face, and back face */

    /*   Create 4D views of receive buffers */
    Vec4D Wrecv(udata->grid->GetRecvView("WEST"), 1, nyl, nzl, dof);
    Vec4D Srecv(udata->grid->GetRecvView("SOUTH"), nxl, 1, nzl, dof);
    Vec4D Brecv(udata->grid->GetRecvView("BACK"), nxl, nyl, 1, dof);

    /*   Perform calculations on each "lower" face */
    Kokkos::parallel_for(
      "AdvectionBoundaryWest", Range3D({0, 0, 0}, {nyl, nzl, dof}),
      KOKKOS_LAMBDA(int j, int k, int l) {
        const int i              = 0;
        const sunrealtype Yijkl  = Yview(i, j, k, l);
        const sunrealtype YSouth = (j > 0) ? Yview(i, j - 1, k, l)
                                           : Srecv(i, 0, k, l);
        const sunrealtype YBack  = (k > 0) ? Yview(i, j, k - 1, l)
                                           : Brecv(i, j, 0, l);
        dYview(i, j, k, l)       = cx * (Yijkl - Wrecv(0, j, k, l)); // d/dx
        dYview(i, j, k, l) += cy * (Yijkl - YSouth);                 // d/dy
        dYview(i, j, k, l) += cz * (Yijkl - YBack);                  // d/dz
      });
    Kokkos::parallel_for(
      "AdvectionBoundarySouth", Range3D({0, 0, 0}, {nxl, nzl, dof}),
      KOKKOS_LAMBDA(int i, int k, int l) {
        const int j             = 0;
        const sunrealtype Yijkl = Yview(i, j, k, l);
        const sunrealtype YWest = (i > 0) ? Yview(i - 1, j, k, l)
                                          : Wrecv(0, j, k, l);
        const sunrealtype YBack = (k > 0) ? Yview(i, j, k - 1, l)
                                          : Brecv(i, j, 0, l);
        dYview(i, j, k, l)      = cx * (Yijkl - YWest);         // d/dx
        dYview(i, j, k, l) += cy * (Yijkl - Srecv(i, 0, k, l)); // d/dy
        dYview(i, j, k, l) += cz * (Yijkl - YBack);             // d/dz
      });
    Kokkos::parallel_for(
      "AdvectionBoundaryBack", Range3D({0, 0, 0}, {nxl, nyl, dof}),
      KOKKOS_LAMBDA(int i, int j, int l) {
        const int k              = 0;
        const sunrealtype Yijkl  = Yview(i, j, k, l);
        const sunrealtype YWest  = (i > 0) ? Yview(i - 1, j, k, l)
                                           : Wrecv(0, j, k, l);
        const sunrealtype YSouth = (j > 0) ? Yview(i, j - 1, k, l)
                                           : Srecv(i, 0, k, l);
        dYview(i, j, k, l)       = cx * (Yijkl - YWest);        // d/dx
        dYview(i, j, k, l) += cy * (Yijkl - YSouth);            // d/dy
        dYview(i, j, k, l) += cz * (Yijkl - Brecv(i, j, 0, l)); // d/dz
      });
  }
  else if (c < 0.0)
  {
    /* Flow moving in the negative x,y,z direction:
       boundaries are east face, north face, and front face */

    /*   Create 4D views of receive buffers */
    Vec4D Erecv(udata->grid->GetRecvView("EAST"), 1, nyl, nzl, dof);
    Vec4D Nrecv(udata->grid->GetRecvView("NORTH"), nxl, 1, nzl, dof);
    Vec4D Frecv(udata->grid->GetRecvView("FRONT"), nxl, nyl, 1, dof);

    /*   Perform calculations on each "upper" face */
    Kokkos::parallel_for(
      "AdvectionBoundaryEast", Range3D({0, 0, 0}, {nyl, nzl, dof}),
      KOKKOS_LAMBDA(int j, int k, int l) {
        const int i              = nxl - 1;
        const sunrealtype Yijkl  = Yview(i, j, k, l);
        const sunrealtype YNorth = (j < nyl - 1) ? Yview(i, j + 1, k, l)
                                                 : Nrecv(i, 0, k, l);
        const sunrealtype YFront = (k < nzl - 1) ? Yview(i, j, k + 1, l)
                                                 : Frecv(i, j, 0, l);
        dYview(i, j, k, l)       = cx * (Erecv(0, j, k, l) - Yijkl); // d/dx
        dYview(i, j, k, l) += cy * (YNorth - Yijkl);                 // d/dy
        dYview(i, j, k, l) += cz * (YFront - Yijkl);                 // d/dz
      });
    Kokkos::parallel_for(
      "AdvectionBoundaryNorth", Range3D({0, 0, 0}, {nxl, nzl, dof}),
      KOKKOS_LAMBDA(int i, int k, int l) {
        const int j              = nyl - 1;
        const sunrealtype Yijkl  = Yview(i, j, k, l);
        const sunrealtype YEast  = (i < nxl - 1) ? Yview(i + 1, j, k, l)
                                                 : Erecv(0, j, k, l);
        const sunrealtype YFront = (k < nzl - 1) ? Yview(i, j, k + 1, l)
                                                 : Frecv(i, j, 0, l);
        dYview(i, j, k, l)       = cx * (YEast - Yijkl);        // d/dx
        dYview(i, j, k, l) += cy * (Nrecv(i, 0, k, l) - Yijkl); // d/dy
        dYview(i, j, k, l) += cz * (YFront - Yijkl);            // d/dz
      });
    Kokkos::parallel_for(
      "AdvectionBoundaryFront", Range3D({0, 0, 0}, {nxl, nyl, dof}),
      KOKKOS_LAMBDA(int i, int j, int l) {
        const int k              = nzl - 1;
        const sunrealtype Yijkl  = Yview(i, j, k, l);
        const sunrealtype YEast  = (i < nxl - 1) ? Yview(i + 1, j, k, l)
                                                 : Erecv(0, j, k, l);
        const sunrealtype YNorth = (j < nyl - 1) ? Yview(i, j + 1, k, l)
                                                 : Nrecv(i, 0, k, l);
        dYview(i, j, k, l)       = cx * (YEast - Yijkl);        // d/dx
        dYview(i, j, k, l) += cy * (YNorth - Yijkl);            // d/dy
        dYview(i, j, k, l) += cz * (Frecv(i, j, 0, l) - Yijkl); // d/dz
      });
  }

  /* return success */
  return (0);
}

/* Compute the reaction term g(t,y). */
static int Reaction(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* set variable shortcuts */
  const sunrealtype A  = udata->A;
  const sunrealtype B  = udata->B;
  const sunrealtype k1 = udata->k1;
  const sunrealtype k2 = udata->k2;
  const sunrealtype k3 = udata->k3;
  const sunrealtype k4 = udata->k4;
  const sunrealtype k5 = udata->k5;
  const sunrealtype k6 = udata->k6;
  const int nxl        = udata->grid->nxl;
  const int nyl        = udata->grid->nyl;
  const int nzl        = udata->grid->nzl;
  const int dof        = udata->grid->dof;

  /* Zero output if not adding reactions to existing RHS */
  if (!udata->add_reactions) { N_VConst(0.0, ydot); }

  /* create 4D views of state and RHS vectors */
  Vec4D Yview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y)), nxl, nyl,
              nzl, dof);
  Vec4D dYview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(ydot)), nxl,
               nyl, nzl, dof);

  /* add reaction terms to RHS */
  Kokkos::parallel_for(
    "ReactionRHS", Range3D({0, 0, 0}, {nxl, nyl, nzl}),
    KOKKOS_LAMBDA(int i, int j, int k) {
      const sunrealtype u = Yview(i, j, k, 0);
      const sunrealtype v = Yview(i, j, k, 1);
      const sunrealtype w = Yview(i, j, k, 2);
      dYview(i, j, k, 0) += k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYview(i, j, k, 1) += k2 * w * u - k3 * u * u * v;
      dYview(i, j, k, 2) += -k2 * w * u + k5 * B - k6 * w;
    });

  /* return success */
  return (0);
}

/* Compute the RHS as h(t,y) = f(t,y) + g(t,y). */
static int AdvectionReaction(sunrealtype t, N_Vector y, N_Vector ydot,
                             void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*)user_data;
  int retval;

  /* NOTE: The order in which Advection and Reaction are called
           is critical here. Advection must be computed first. */
  retval = Advection(t, y, ydot, user_data);
  if (check_retval((void*)&retval, "Advection", 1, udata->myid))
  {
    return (-1);
  }

  retval = Reaction(t, y, ydot, user_data);
  if (check_retval((void*)&retval, "Reaction", 1, udata->myid)) { return (-1); }

  /* return success */
  return (0);
}

/* Compute the residual F(t,y,y') = ydot - h(t,y) = 0. */
static int AdvectionReactionResidual(sunrealtype t, N_Vector y, N_Vector ydot,
                                     N_Vector F, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*)user_data;
  int retval;

  /* NOTE: The order in which Advection and Reaction are called
           is critical here. Advection must be computed first. */
  retval = Advection(t, y, F, user_data); /* F = -c y_x */
  if (check_retval((void*)&retval, "Advection", 1, udata->myid))
  {
    return (-1);
  }

  retval = Reaction(t, y, F, user_data); /* F = -c y_x + g(t,y) */
  if (check_retval((void*)&retval, "Reaction", 1, udata->myid)) { return (-1); }

  /* F = ydot - h(t,y) = ydot + c y_x - g(t,y) */
  N_VLinearSum(1.0, ydot, -1.0, F, F);

  /* return success */
  return (0);
}

/* --------------------------------------------------------------
 * Linear system and Jacobian functions
 * --------------------------------------------------------------*/

/* Solve the linear systems Ax = b where A = I - gamma*dg/dy.
   When using a fully implicit method, we are approximating
   dh/dy as dg/dy. */
static int SolveReactionLinSys(N_Vector y, N_Vector x, N_Vector b,
                               const sunrealtype gamma, UserData* udata)
{
  /* set variable shortcuts */
  const int dof        = udata->grid->dof;
  const int nxl        = udata->grid->nxl;
  const int nyl        = udata->grid->nyl;
  const int nzl        = udata->grid->nzl;
  const sunrealtype k2 = udata->k2;
  const sunrealtype k3 = udata->k3;
  const sunrealtype k4 = udata->k4;
  const sunrealtype k6 = udata->k6;

  /* create 4D views of state, RHS and solution vectors */
  Vec4D Yview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y)), nxl, nyl,
              nzl, dof);
  Vec4D Bview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(b)), nxl, nyl,
              nzl, dof);
  Vec4D Xview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(x)), nxl, nyl,
              nzl, dof);

  /* solve reaction linear system */
  Kokkos::parallel_for(
    "SolveReactionLinSys", Range3D({0, 0, 0}, {nxl, nyl, nzl}),
    KOKKOS_LAMBDA(int i, int j, int k) {
      /* shortcuts to u, v, w for the block */
      const sunrealtype u = Yview(i, j, k, 0);
      const sunrealtype v = Yview(i, j, k, 1);
      const sunrealtype w = Yview(i, j, k, 2);

      //
      // compute A = I - gamma*(dg/dy)
      //

      /* 1st row: u, v, w */
      const sunrealtype A0 = 1. - gamma * (-k2 * w + 2.0 * k3 * u * v - k4);
      const sunrealtype A1 = -gamma * (k3 * u * u);
      const sunrealtype A2 = -gamma * (-k2 * u);

      /* 2nd row: u, v, w */
      const sunrealtype A3 = -gamma * (k2 * w - 2.0 * k3 * u * v);
      const sunrealtype A4 = 1. - gamma * (-k3 * u * u);
      const sunrealtype A5 = -gamma * (k2 * u);

      /* 3rd row: u, v, w */
      const sunrealtype A6 = -gamma * (-k2 * w);
      const sunrealtype A7 = 0.0;
      const sunrealtype A8 = 1. - gamma * (-k2 * u - k6);

      //
      // compute x = A^{-1}*b
      //

      const sunrealtype scratch_0  = A4 * A8;
      const sunrealtype scratch_1  = A1 * A5;
      const sunrealtype scratch_2  = A2 * A7;
      const sunrealtype scratch_3  = A5 * A7;
      const sunrealtype scratch_4  = A1 * A8;
      const sunrealtype scratch_5  = A2 * A4;
      const sunrealtype scratch_6  = 1.0 / (A0 * scratch_0 - A0 * scratch_3 +
                                           A3 * scratch_2 - A3 * scratch_4 +
                                           A6 * scratch_1 - A6 * scratch_5);
      const sunrealtype scratch_7  = A2 * A3;
      const sunrealtype scratch_8  = A6 * Bview(i, j, k, 0);
      const sunrealtype scratch_9  = A2 * A6;
      const sunrealtype scratch_10 = A3 * Bview(i, j, k, 0);
      const sunrealtype scratch_11 = 1.0 / A0;
      const sunrealtype scratch_12 = A1 * scratch_11;
      const sunrealtype scratch_13 = (-A6 * scratch_12 + A7) /
                                     (-A3 * scratch_12 + A4);

      Xview(i, j, k, 0) = scratch_6 *
                          (Bview(i, j, k, 0) * (scratch_0 - scratch_3) +
                           Bview(i, j, k, 1) * (scratch_2 - scratch_4) +
                           Bview(i, j, k, 2) * (scratch_1 - scratch_5));
      Xview(i, j, k, 1) = scratch_6 * (Bview(i, j, k, 2) * (scratch_7 - A0 * A5) +
                                       Bview(i, j, k, 1) * (A0 * A8 - scratch_9) +
                                       A5 * scratch_8 - A8 * scratch_10);
      Xview(i, j, k,
            2) = (-Bview(i, j, k, 2) + scratch_11 * scratch_8 +
                  scratch_13 * (Bview(i, j, k, 1) - scratch_10 * scratch_11)) /
                 (-A8 + scratch_11 * scratch_9 +
                  scratch_13 * (A5 - scratch_11 * scratch_7));
    });

  return (0);
}

/* Solve the linear systems Ax = b where A = -dg/dy + gamma.
   We are approximating dh/dy as dg/dy. */
static int SolveReactionLinSysRes(N_Vector y, N_Vector x, N_Vector b,
                                  const sunrealtype gamma, UserData* udata)
{
  /* set variable shortcuts */
  const int dof        = udata->grid->dof;
  const int nxl        = udata->grid->nxl;
  const int nyl        = udata->grid->nyl;
  const int nzl        = udata->grid->nzl;
  const sunrealtype k2 = udata->k2;
  const sunrealtype k3 = udata->k3;
  const sunrealtype k4 = udata->k4;
  const sunrealtype k6 = udata->k6;

  /* create 4D views of state, RHS and solution vectors */
  Vec4D Yview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y)), nxl, nyl,
              nzl, dof);
  Vec4D Bview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(b)), nxl, nyl,
              nzl, dof);
  Vec4D Xview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(x)), nxl, nyl,
              nzl, dof);

  /* solve reaction linear system */
  Kokkos::parallel_for(
    "SolveReactionLinSys", Range3D({0, 0, 0}, {nxl, nyl, nzl}),
    KOKKOS_LAMBDA(int i, int j, int k) {
      /* shortcuts to u, v, w for the block */
      const sunrealtype u = Yview(i, j, k, 0);
      const sunrealtype v = Yview(i, j, k, 1);
      const sunrealtype w = Yview(i, j, k, 2);

      //
      // compute A = -dg/dy + gamma*diag(df/dydot)
      // where diag(df/dydot) is approximated as
      // diag([udot, vdot, wdot])
      //

      /* 1st row: u, v, w */
      const sunrealtype A0 = -(-k2 * w + 2.0 * k3 * u * v - k4) + gamma;
      const sunrealtype A1 = -(k3 * u * u);
      const sunrealtype A2 = -(-k2 * u);

      /* 2nd row: u, v, w */
      const sunrealtype A3 = -(k2 * w - 2.0 * k3 * u * v);
      const sunrealtype A4 = -(-k3 * u * u) + gamma;
      const sunrealtype A5 = -(k2 * u);

      /* 3rd row: u, v, w */
      const sunrealtype A6 = -(-k2 * w);
      const sunrealtype A7 = 0.0;
      const sunrealtype A8 = -(-k2 * u - k6) + gamma;

      //
      // compute x = A^{-1}*b
      //

      const sunrealtype scratch_0  = A4 * A8;
      const sunrealtype scratch_1  = A1 * A5;
      const sunrealtype scratch_2  = A2 * A7;
      const sunrealtype scratch_3  = A5 * A7;
      const sunrealtype scratch_4  = A1 * A8;
      const sunrealtype scratch_5  = A2 * A4;
      const sunrealtype scratch_6  = 1.0 / (A0 * scratch_0 - A0 * scratch_3 +
                                           A3 * scratch_2 - A3 * scratch_4 +
                                           A6 * scratch_1 - A6 * scratch_5);
      const sunrealtype scratch_7  = A2 * A3;
      const sunrealtype scratch_8  = A6 * Bview(i, j, k, 0);
      const sunrealtype scratch_9  = A2 * A6;
      const sunrealtype scratch_10 = A3 * Bview(i, j, k, 0);
      const sunrealtype scratch_11 = 1.0 / A0;
      const sunrealtype scratch_12 = A1 * scratch_11;
      const sunrealtype scratch_13 = (-A6 * scratch_12 + A7) /
                                     (-A3 * scratch_12 + A4);

      Xview(i, j, k, 0) = scratch_6 *
                          (Bview(i, j, k, 0) * (scratch_0 - scratch_3) +
                           Bview(i, j, k, 1) * (scratch_2 - scratch_4) +
                           Bview(i, j, k, 2) * (scratch_1 - scratch_5));
      Xview(i, j, k, 1) = scratch_6 * (Bview(i, j, k, 2) * (scratch_7 - A0 * A5) +
                                       Bview(i, j, k, 1) * (A0 * A8 - scratch_9) +
                                       A5 * scratch_8 - A8 * scratch_10);
      Xview(i, j, k,
            2) = (-Bview(i, j, k, 2) + scratch_11 * scratch_8 +
                  scratch_13 * (Bview(i, j, k, 1) - scratch_10 * scratch_11)) /
                 (-A8 + scratch_11 * scratch_9 +
                  scratch_13 * (A5 - scratch_11 * scratch_7));
    });

  return (0);
}

/* --------------------------------------------------------------
 * Preconditioner functions
 * --------------------------------------------------------------*/

/* Solves Pz = r where P = I - gamma * dg/dy */
static int PSolve(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector r,
                  N_Vector z, sunrealtype gamma, sunrealtype delta, int lr,
                  void* user_data)
{
  /* local variables */
  UserData* udata = (UserData*)user_data;
  int retval;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* solve the task-local linear system Pz = r */
  retval = SolveReactionLinSys(y, z, r, gamma, udata);

  return (retval);
}

/* Solves Pz = r where P = -dg/dy + gamma */
static int PSolveRes(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector F,
                     N_Vector r, N_Vector z, sunrealtype cj, sunrealtype delta,
                     void* user_data)
{
  /* local variables */
  UserData* udata = (UserData*)user_data;
  int retval;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* solve the task-local linear system Pz = r */
  retval = SolveReactionLinSysRes(y, z, r, cj, udata);

  return (retval);
}

#endif
