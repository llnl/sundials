/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * -----------------------------------------------------------------------------
 * Jacobi preconditiner for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

#if defined(BENCHMARK_ODE)

// Preconditioner setup routine
int PSetup(sunrealtype t, N_Vector u, N_Vector f, sunbooleantype jok,
           sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Constants for computing diffusion
  sunrealtype cx = udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = -TWO * (cx + cy);

  // Set all entries of d to the inverse diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  sunrealtype c = ONE / (ONE - gamma * cc);
  N_VConst(c, udata->diag);

  // Return success
  return 0;
}

// Preconditioner solve routine for Pz = r
int PSolve(sunrealtype t, N_Vector u, N_Vector f, N_Vector r, N_Vector z,
           sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  // Access user_data structure
  UserData* udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Perform Jacobi iteration
  N_VProd(udata->diag, r, z);

  // Return success
  return 0;
}

#elif defined(BENCHMARK_DAE)

// Preconditioner setup and solve functions
int PSetup(sunrealtype t, N_Vector u, N_Vector up, N_Vector res, sunrealtype cj,
           void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Constants for computing diffusion
  sunrealtype cx = udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = -TWO * (cx + cy);

  // Set all entries of d to the inverse diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  sunrealtype c = ONE / (cj - cc);
  N_VConst(c, udata->diag);

  // Return success
  return 0;
}

int PSolve(sunrealtype t, N_Vector u, N_Vector up, N_Vector res, N_Vector r,
           N_Vector z, sunrealtype cj, sunrealtype delta, void* user_data)
{
  // Access user_data structure
  UserData* udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Perform Jacobi iteration
  N_VProd(udata->diag, r, z);

  // Return success
  return 0;
}

#else
#error "Missing ODE/DAE preprocessor directive"
#endif
