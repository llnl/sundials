/*
 * -----------------------------------------------------------------
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
 * This is the header file for a simple C-language math library. The
 * routines listed here work with the type sunrealtype as defined in
 * the header file sundials_types.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALSMATH_H
#define _SUNDIALSMATH_H

#include <math.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Macros
 * -----------------------------------------------------------------
 * SUNMIN(A,B) returns the minimum of A and B
 *
 * SUNMAX(A,B) returns the maximum of A and B
 *
 * SUNSQR(A) returns A^2
 *
 * SUNRsqrt calls the appropriate version of sqrt
 *
 * SUNRabs calls the appropriate version of abs
 *
 * SUNRexp calls the appropriate version of exp
 *
 * SUNRceil calls the appropriate version of ceil
 *
 * SUNRround calls the appropriate version of round
 * -----------------------------------------------------------------
 */

#ifndef SUNMIN
#define SUNMIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUNMAX
#define SUNMAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUNSQR
#define SUNSQR(A) ((A) * (A))
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsqrt
 * -----------------------------------------------------------------
 * Usage : sunrealtype sqrt_x;
 *         sqrt_x = SUNRsqrt(x);
 * -----------------------------------------------------------------
 * SUNRsqrt(x) returns the square root of x. If x < ZERO, then
 * SUNRsqrt returns ZERO.
 * -----------------------------------------------------------------
 */

#ifndef SUNRsqrt
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrt((x))))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrtf((x))))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrtl((x))))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRabs
 * -----------------------------------------------------------------
 * Usage : sunrealtype abs_x;
 *         abs_x = SUNRabs(x);
 * -----------------------------------------------------------------
 * SUNRabs(x) returns the absolute value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRabs
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRabs(x) (fabs((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRabs(x) (fabsf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRabs(x) (fabsl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRexp
 * -----------------------------------------------------------------
 * Usage : sunrealtype exp_x;
 *         exp_x = SUNRexp(x);
 * -----------------------------------------------------------------
 * SUNRexp(x) returns e^x (base-e exponential function).
 * -----------------------------------------------------------------
 */

#ifndef SUNRexp
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRexp(x) (exp((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRexp(x) (expf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRexp(x) (expl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRceil
 * -----------------------------------------------------------------
 * Usage : sunrealtype ceil_x;
 *         ceil_x = SUNRceil(x);
 * -----------------------------------------------------------------
 * SUNRceil(x) returns the smallest integer value not less than x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRceil
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRceil(x) (ceil((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRceil(x) (ceilf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRceil(x) (ceill((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRcopysign
 * -----------------------------------------------------------------
 * Usage : sunrealtype z;
 *         z = SUNRcopysign(x, y);
 * -----------------------------------------------------------------
 * SUNRcopysign(x, y) returns x with the sign of y.
 * -----------------------------------------------------------------
 */

#ifndef SUNRcopysign
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRcopysign(x, y) (copysign((x), (y)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRcopysign(x, y) (copysignf((x), (y)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRcopysign(x, y) (copysignl((x), (y)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsamesign
 * -----------------------------------------------------------------
 * Usage : sunrealtype z;
 *         z = SUNRsamesign(x, y);
 * -----------------------------------------------------------------
 * SUNRsamesign(x, y) returns true if x and y share the same sign,
 * false otherwise
 * -----------------------------------------------------------------
 */

#ifndef SUNRsamesign
#define SUNRsamesign(x, y) (signbit((x)) == signbit((y)))
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRdifferentsign
 * -----------------------------------------------------------------
 * Usage : sunrealtype z;
 *         z = SUNRdifferentsign(x, y);
 * -----------------------------------------------------------------
 * SUNRdifferentsign(x) returns true if x and y have different
 * signs, false otherwise
 * -----------------------------------------------------------------
 */

#ifndef SUNRdifferentsign
#define SUNRdifferentsign(x, y) (!SUNRsamesign((x), (y)))
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerR
 * -----------------------------------------------------------------
 * Usage : sunrealtype base, exponent, ans;
 *         ans = SUNRpowerR(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerR returns the value of base^exponent, where both base and
 * exponent are of type sunrealtype.
 * -----------------------------------------------------------------
 */
#ifndef SUNRpowerR
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRpowerR(base, exponent) (pow(base, exponent))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRpowerR(base, exponent) (powf(base, exponent))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRpowerR(base, exponent) (powl(base, exponent))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRround
 * -----------------------------------------------------------------
 * Usage : sunrealtype round_x;
 *         round_x = SUNRround(x);
 * -----------------------------------------------------------------
 * SUNRround(x) returns the smallest integer value not less than x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRround
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRround(x) (round((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRround(x) (roundf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRround(x) (roundl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNIpowerI
 * -----------------------------------------------------------------
 * Usage : int exponent, base, ans;
 *         ans = SUNIpowerI(base,exponent);
 * -----------------------------------------------------------------
 * SUNIpowerI returns the value of base^exponent, where base and
 * exponent are of type int and exponent is nonnegative.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SUNIpowerI(int base, int exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerI
 * -----------------------------------------------------------------
 * Usage : int exponent;
 *         sunrealtype base, ans;
 *         ans = SUNRpowerI(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerI returns the value of base^exponent, where base is of type
 * sunrealtype and exponent is of type int.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNRpowerI(sunrealtype base, int exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompare
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype a, b;
 *         isNotEqual = SUNRCompare(a, b);
 * -----------------------------------------------------------------
 * SUNRCompareTol returns 0 if the relative difference of a and b is
 * less than or equal to 10*machine epsilon. If the relative
 * difference is greater than 10*machine epsilon, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompareTol
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype a, b, tol;
 *         isNotEqual = SUNRCompareTol(a, b, tol);
 * -----------------------------------------------------------------
 * SUNRCompareTol returns 0 if the relative difference of a and b is
 * less than or equal to the provided tolerance. If the relative
 * difference is greater than the tolerance, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b,
                                              sunrealtype tol);

/*
 * -----------------------------------------------------------------
 * Function : SUNStrToReal
 * -----------------------------------------------------------------
 * Usage : sunrealtype a = SUNStrToReal(const char* str)
 * -----------------------------------------------------------------
 * SUNStrToReal parses str into the sunrealtype variable. Uses standard
 * strtod variants when they are available.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNStrToReal(const char* str);

#ifdef __cplusplus
}
#endif

#endif
