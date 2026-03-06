/* -----------------------------------------------------------------------------
 * Programmer(s): Steven Roberts @ LLNL
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
 * Unit test that checks splitting coefficients constructors and loading
 * functions
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_splittingstep.h"
#include "sundials/sundials_math.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)
#define HALF SUN_RCONST(0.5)
#define GAMMA         \
  (SUN_RCONST(1.0) /  \
   (SUN_RCONST(4.0) - \
    SUNRpowerR(SUN_RCONST(4.0), SUN_RCONST(1.0) / SUN_RCONST(3.0))))

static int check_coefficients(const char* name,
                              SplittingStepCoefficients coefficients,
                              int expected_sequential_methods,
                              int expected_stages, int expected_partitions,
                              int expected_order,
                              const sunrealtype* expected_alpha,
                              const sunrealtype* expected_beta)
{
  printf("Testing %s\n", name);

  if (coefficients == NULL)
  {
    fprintf(stderr, "Expected non-NULL coefficients\n");
  }

  if (coefficients->sequential_methods != expected_sequential_methods)
  {
    fprintf(stderr, "Expected %d sequential methods, but there are %d\n",
            expected_sequential_methods, coefficients->sequential_methods);
    return 1;
  }

  if (coefficients->stages != expected_stages)
  {
    fprintf(stderr, "Expected %d stages, but there are %d\n", expected_stages,
            coefficients->stages);
    return 1;
  }

  if (coefficients->partitions != expected_partitions)
  {
    fprintf(stderr, "Expected %d partitions, but there are %d\n",
            expected_partitions, coefficients->partitions);
    return 1;
  }

  if (coefficients->order != expected_order)
  {
    fprintf(stderr, "Expected order %d, but is %d\n", expected_order,
            coefficients->order);
    return 1;
  }

  int retval = 0;

  for (int i = 0; i < expected_sequential_methods; i++)
  {
    if (SUNRCompare(coefficients->alpha[i], expected_alpha[i]))
    {
      fprintf(stderr, "alpha[%d] incorrect\n", i);
      retval++;
    }
  }

  for (int i = 0; i < expected_sequential_methods; i++)
  {
    for (int j = 0; j <= expected_stages; j++)
    {
      for (int k = 0; k < expected_partitions; k++)
      {
        if (SUNRCompare(coefficients->beta[i][j][k],
                        expected_beta[(i * (expected_stages + 1) + j) * expected_partitions +
                                      k]))
        {
          fprintf(stderr, "beta[%d][%d][%d] incorrect\n", i, j, k);
          retval++;
        }
      }
    }
  }

  if (retval > 0) { SplittingStepCoefficients_Write(coefficients, stderr); }

  return retval;
}

/* Main program */
int main(int argc, char* argv[])
{
  int retval = 0;

  printf("Testing Splitting Coefficients\n");

  SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Create(1, 2, 2, 1, NULL, NULL);
  if (coefficients != NULL)
  {
    fprintf(stderr, "Coefficients created with NULL coefficients\n");
    retval++;
  }

  sunrealtype alpha_lie_trotter[] = {ONE};
  sunrealtype beta_lie_trotter[]  = {ZERO, ZERO, ONE, ONE};

  coefficients = SplittingStepCoefficients_Create(0, 0, 0, 1, alpha_lie_trotter,
                                                  beta_lie_trotter);
  if (coefficients != NULL)
  {
    fprintf(stderr, "Coefficients created with invalid sizes\n");
    retval++;
  }

  coefficients = SplittingStepCoefficients_Create(1, 1, 2, 1, alpha_lie_trotter,
                                                  beta_lie_trotter);
  retval += check_coefficients("Lie-Trotter (manually created)", coefficients,
                               1, 1, 2, 1, alpha_lie_trotter, beta_lie_trotter);

  SplittingStepCoefficients coefficients_copy =
    SplittingStepCoefficients_Copy(coefficients);
  retval += check_coefficients("Lie-Trotter (copy)", coefficients_copy, 1, 1, 2,
                               1, alpha_lie_trotter, beta_lie_trotter);
  SplittingStepCoefficients_Destroy(&coefficients);
  SplittingStepCoefficients_Destroy(&coefficients_copy);

  coefficients = SplittingStepCoefficients_LoadCoefficients(
    ARKODE_SPLITTING_LIE_TROTTER_1_1_2);
  retval += check_coefficients("Lie-Trotter (load by enum)", coefficients, 1, 1,
                               2, 1, alpha_lie_trotter, beta_lie_trotter);
  SplittingStepCoefficients_Destroy(&coefficients);

  coefficients = SplittingStepCoefficients_LoadCoefficientsByName(
    "ARKODE_SPLITTING_LIE_TROTTER_1_1_2");
  retval += check_coefficients("Lie-Trotter (load by name)", coefficients, 1, 1,
                               2, 1, alpha_lie_trotter, beta_lie_trotter);
  SplittingStepCoefficients_Destroy(&coefficients);

  coefficients = SplittingStepCoefficients_LieTrotter(2);
  retval += check_coefficients("Lie-Trotter (constructor)", coefficients, 1, 1,
                               2, 1, alpha_lie_trotter, beta_lie_trotter);
  SplittingStepCoefficients_Destroy(&coefficients);

  sunrealtype alpha_strang[] = {ONE};
  sunrealtype beta_strang[]  = {ZERO, ZERO, ZERO, HALF, HALF, ONE,
                                HALF, ONE,  ONE,  ONE,  ONE,  ONE};

  coefficients = SplittingStepCoefficients_Strang(3);
  retval += check_coefficients("Strang", coefficients, 1, 3, 3, 2, alpha_strang,
                               beta_strang);
  SplittingStepCoefficients_Destroy(&coefficients);

  sunrealtype alpha_parallel[] = {ONE, ONE, -ONE};
  sunrealtype beta_parallel[]  = {ZERO, ZERO, ONE,  ZERO, ZERO, ZERO,
                                  ZERO, ONE,  ZERO, ZERO, ZERO, ZERO};

  coefficients = SplittingStepCoefficients_Parallel(2);
  retval += check_coefficients("Parallel", coefficients, 3, 1, 2, 1,
                               alpha_parallel, beta_parallel);
  SplittingStepCoefficients_Destroy(&coefficients);

  sunrealtype alpha_symmetric_parallel[] = {HALF, HALF, -ONE};
  sunrealtype beta_symmetric_parallel[]  = {ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,
                                            ZERO, ZERO, ZERO, ONE,  ONE,  ONE,
                                            ZERO, ZERO, ZERO, ZERO, ZERO, ONE,
                                            ZERO, ONE,  ONE,  ONE,  ONE,  ONE};

  coefficients = SplittingStepCoefficients_SymmetricParallel(3);
  retval += check_coefficients("Symmetric Parallel", coefficients, 2, 3, 3, 2,
                               alpha_symmetric_parallel, beta_symmetric_parallel);
  SplittingStepCoefficients_Destroy(&coefficients);

  sunrealtype alpha_suzuki_fractal[] = {ONE};
  sunrealtype beta_suzuki_fractal[]  = {ZERO,
                                        ZERO,
                                        HALF * GAMMA,
                                        GAMMA,
                                        (ONE + HALF) * GAMMA,
                                        TWO * GAMMA,
                                        HALF,
                                        ONE - TWO * GAMMA,
                                        ONE - (ONE + HALF) * GAMMA,
                                        ONE - GAMMA,
                                        ONE - HALF * GAMMA,
                                        ONE,
                                        ONE,
                                        ONE};

  coefficients = SplittingStepCoefficients_SuzukiFractal(2, 4);
  retval += check_coefficients("Suzuki Fractal", coefficients, 1, 6, 2, 4,
                               alpha_suzuki_fractal, beta_suzuki_fractal);
  SplittingStepCoefficients_Destroy(&coefficients);

  printf("%d test failures\n", retval);

  return retval;
}
