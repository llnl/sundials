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
 * Unit test for reusing a KINSOL instance based on bug report in the MFEM repo
 * https://github.com/mfem/mfem/issues/5004
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kinsol/kinsol.h"
#include "nvector/nvector_serial.h"

/* Fixed point function g(x) = 0.5x^2 + 0.25
 * x = 0.5 x^2 + 0.25 -> roots ~ 0.292893 and 1.707106.
 * At x* ~ 0.292893, g'(x*) = x* < 1 => contraction. */
static int gfun(N_Vector x, N_Vector g, void* user_data)
{
  sunrealtype* x_data = N_VGetArrayPointer(x);
  if (!x_data) { return -1; }
  sunrealtype* g_data = N_VGetArrayPointer(g);
  if (!g_data) { return -1; }

  const double xv = x_data[0];
  g_data[0]       = SUN_RCONST(0.5) * xv * xv + SUN_RCONST(0.25);

  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  int retval        = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval)
  {
    printf("ERROR: SUNContext_Create returned %i", retval);
    return 1;
  }

  // Initial guess/solution vector
  N_Vector x = N_VNew_Serial(1, sunctx);
  if (!x)
  {
    printf("ERROR: N_VNew_Serial returned NULL");
    return 1;
  }

  // Scaling vector
  N_Vector scale = N_VClone(x);
  if (!scale)
  {
    printf("ERROR: N_VClone returned NULL");
    return 1;
  }
  N_VConst(1.0, scale);

  // Create KINSOL
  void* kmem = KINCreate(sunctx);
  if (!kmem)
  {
    printf("ERROR: KINCreate returned NULL");
    return 1;
  }

  retval = KINInit(kmem, gfun, x);
  if (retval)
  {
    printf("ERROR: KINInit returned %i", retval);
    return 1;
  }

  // Set options
  int maa = 0;
  if (argc > 1) { maa = atoi(argv[1]); }
  printf("Anderson acceleration depth = %i\n", maa);
  retval = KINSetMAA(kmem, maa);
  if (retval)
  {
    printf("ERROR: KINSetMAA returned %i", retval);
    return 1;
  }

  // -------
  // Solve 1
  // -------

  // Set the initial guess
  sunrealtype* x_data = N_VGetArrayPointer(x);
  if (!x_data)
  {
    printf("ERROR: N_VGetArrayPointer returned NULL");
    return 1;
  }
  x_data[0] = SUN_RCONST(0.10);

  // Solve the system
  retval = KINSol(kmem, x, KIN_FP, scale, scale);
  if (retval)
  {
    printf("ERROR: KINSol returned %i", retval);
    return 1;
  }

  printf("\nSolution 1:\n");
  N_VPrint(x);
  sunrealtype sol_1 = x_data[0];

  // Get solver stats
  long int nni_1;
  retval = KINGetNumNonlinSolvIters(kmem, &nni_1);
  if (retval)
  {
    printf("ERROR: KINGetNumNonlinSolvIters returned %i", retval);
    return 1;
  }

  long int nfe_1;
  retval = KINGetNumFuncEvals(kmem, &nfe_1);
  if (retval)
  {
    printf("ERROR: KINGetNumFuncEvals returned %i", retval);
    return 1;
  }

  printf("\nFinal Statistics:\n");
  printf("Number of nonlinear iterations: %6ld\n", nni_1);
  printf("Number of function evaluations: %6ld\n", nfe_1);

  // -------
  // Solve 2
  // -------

  // Reset the initial guess
  x_data[0] = SUN_RCONST(0.10);

  // Solve the system
  retval = KINSol(kmem, x, KIN_FP, scale, scale);
  if (retval)
  {
    printf("ERROR: KINSol returned %i", retval);
    return 1;
  }

  printf("\nSolution 2:\n");
  N_VPrint(x);
  sunrealtype sol_2 = x_data[0];

  // Get solver stats
  long int nni_2;
  retval = KINGetNumNonlinSolvIters(kmem, &nni_2);
  if (retval)
  {
    printf("ERROR: KINGetNumNonlinSolvIters returned %i", retval);
    return 1;
  }

  long int nfe_2;
  retval = KINGetNumFuncEvals(kmem, &nfe_2);
  if (retval)
  {
    printf("ERROR: KINGetNumFuncEvals returned %i", retval);
    return 1;
  }

  printf("\nFinal Statistics:\n");
  printf("Number of nonlinear iterations: %6ld\n", nni_2);
  printf("Number of function evaluations: %6ld\n", nfe_2);

  // --------------
  // Compare solves
  // --------------

  // Solutions should be identical
  if (sol_1 != sol_2)
  {
    printf("\nERROR: Solutions differ!\n");
    printf("Solution 1: " SUN_FORMAT_G "\n", sol_1);
    printf("Solution 2: " SUN_FORMAT_G "\n", sol_2);
    return 1;
  }

  // Iterations should be identical
  if (nni_1 != nni_2)
  {
    printf("\nERROR: Number of iterations differ!\n");
    printf("Iterations 1: %ld\n", nni_1);
    printf("Iterations 2: %ld\n", nni_2);
    return 1;
  }

  // Evaluations should be identical
  if (nfe_1 != nfe_2)
  {
    printf("\nERROR: Number of function evaluations differ!\n");
    printf("Evaluations 1: %ld\n", nfe_1);
    printf("Evaluations 2: %ld\n", nfe_2);
    return 1;
  }

  // -----------
  // Free memory
  // -----------

  N_VDestroy(x);
  N_VDestroy(scale);
  KINFree(&kmem);
  SUNContext_Free(&sunctx);

  return 0;
}
