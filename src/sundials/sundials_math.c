/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the implementation file for a simple C-language math
 * library.
 * -----------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

int SUNIpowerI(int base, int exponent)
{
  int i;
  int prod = 1;

  for (i = 1; i <= exponent; i++) { prod *= base; }
  return (prod);
}

sunrealtype SUNRpowerI(sunrealtype base, int exponent)
{
  int i, expt;
  sunrealtype prod;

  prod = SUN_RCONST(1.0);
  expt = abs(exponent);
  for (i = 1; i <= expt; i++) { prod *= base; }
  if (exponent < 0) { prod = SUN_RCONST(1.0) / prod; }
  return (prod);
}

sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b)
{
  return (SUNRCompareTol(a, b, 10 * SUN_UNIT_ROUNDOFF));
}

sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol)
{
  sunrealtype diff;
  sunrealtype norm;

  /* If a and b are exactly equal.
   * This also covers the case where a and b are both inf under IEEE 754.
   */
  if (a == b) { return (SUNFALSE); }

  diff = SUNRabs(a - b);
  norm = SUNMIN(SUNRabs(a + b), SUN_BIG_REAL);

  /* When |a + b| is very small (less than 10*SUN_UNIT_ROUNDOFF) or zero, we use
   * an absolute difference:
   *    |a - b| >= 10*SUN_UNIT_ROUNDOFF
   * Otherwise we use a relative difference:
   *    |a - b| < tol * |a + b|
   * The choice to use |a + b| over max(a, b) is arbitrary, as is the choice to
   * use 10*SUN_UNIT_ROUNDOFF.
   * 
   * In order to handle NANs correctly without explicit checks of isnan or
   * isunordered (which throw warnings for some compilers and flags), we use
   * !isless. The seemingly equivalent >= can have undefined behavior for NANs.
   */
  return !isless(diff, SUNMAX(10 * SUN_UNIT_ROUNDOFF, tol * norm));
}

sunrealtype SUNStrToReal(const char* str)
{
  char* end;
#if defined(SUNDIALS_EXTENDED_PRECISION)
  return strtold(str, &end);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  return strtod(str, &end);
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return strtof(str, &end);
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
}
