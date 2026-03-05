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
#ifdef SUNDIALS_FLOAT128_PRECISION
/* This defines an output stream operator for the `__float128` type.*/
#include <iostream>
#include <iomanip>
#include <quadmath.h>
#include <cstdio>

static std::ostream& operator<<(std::ostream& os, __float128 value)
{
  // Get current stream formatting state
  const int width = os.width();     // Width set by std::setw
  const int precision = os.precision(); // Precision set by std::setprecision
  const std::ios_base::fmtflags flags = os.flags(); // Format flags (e.g., scientific notation)

  // Determine format specifier based on stream flags (e/f/g)
  char format_specifier = 'g';
  if (flags & std::ios_base::scientific)
  {
    format_specifier = 'e';
  }
  else if (flags & std::ios_base::fixed)
  {
    format_specifier = 'f';
  }

  // Dynamically generate format string (e.g., "%20.15Qe")
  char format_buffer[64];
  std::snprintf(
    format_buffer, sizeof(format_buffer),
    "%%%d.%dQ%c",  // Format template: %[width].[precision]Q[e/f/g]
    width,         // Width from setw
    precision,     // Precision from setprecision
    format_specifier
  );

  // Format __float128 to string
  char value_buffer[128];
  int n = quadmath_snprintf(
    value_buffer, sizeof(value_buffer),
    format_buffer, // Dynamically generated format (e.g., "%20.15Qe")
    value
  );

  // Write to output stream
  if (n >= 0 && n < sizeof(value_buffer))
  {
    os << value_buffer;
  }
  else
  {
    os << "[FORMAT ERROR]";
  }

  // Reset stream width (setw has one-time effect)
  os.width(0);

  return os;
}

#endif

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
 * SUNRisnan calls the appropriate version of isnan
 *
 * SUNRexp calls the appropriate version of exp
 *
 * SUNRlog calls the appropriate version of log
 *
 * SUNRceil calls the appropriate version of ceil
 *
 * SUNRcopysign calls the appropriate version copysign
 *
 * SUNRpowerR calls the appropriate version pow
 *
 * SUNRround calls the appropriate version of round
 *
 * SUNRsin calls the appropriate version of sin
 *
 * SUNRcos calls the appropriate version of cos
 *
 * SUNRasin calls the appropriate version of asin
 *
 * SUNRacos calls the appropriate version of acos
 *
 * SUNRsinh calls the appropriate version of sinh
 *
 * SUNRcosh calls the appropriate version of cosh

 * SUNRatan calls the appropriate version of atan
 *
 * SUNRpower calls the appropriate version of power
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRsqrt(x) ((x) <= SUN_RCONST(0.0) ? (SUN_RCONST(0.0)) : (sqrtq((x))))
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRabs(x) (fabsq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRisnan
 * -----------------------------------------------------------------
 * Usage : sunrealtype isnan_x;
 *         isnan_x = SUNRisnan(x);
 * -----------------------------------------------------------------
 * SUNRisnan(x) returns isnan_x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRisnan
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRisnan(x) (isnan((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRisnan(x) (isnanf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRisnan(x) (isnanl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRisnan(x) (isnanq((x)))
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRexp(x) (expq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRlog
 * -----------------------------------------------------------------
 * Usage : sunrealtype log_x;
 *         log_x = SUNRlog(x);
 * -----------------------------------------------------------------
 * SUNRlog(x) returns log_x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRlog
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRlog(x) (log((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRlog(x) (logf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRlog(x) (logl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRlog(x) (logq((x)))
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRceil(x) (ceilq((x)))
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRcopysign(x, y) (copysignq((x), (y)))
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRpowerR(base, exponent) (powq(base, exponent))
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
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRround(x) (roundq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsin
 * -----------------------------------------------------------------
 * Usage : sunrealtype sin_x;
 *         sin_x = SUNRsin(x);
 * -----------------------------------------------------------------
 * SUNRsin(x) returns the sin value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRsin
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsin(x) (sin((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsin(x) (sinf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsin(x) (sinl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRsin(x) (sinq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRcos
 * -----------------------------------------------------------------
 * Usage : sunrealtype cos_x;
 *         cos_x = SUNRcos(x);
 * -----------------------------------------------------------------
 * SUNRcos(x) returns the cos value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRcos
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRcos(x) (cos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRcos(x) (cosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRcos(x) (cosl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRcos(x) (cosq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRasin
 * -----------------------------------------------------------------
 * Usage : sunrealtype cos_x;
 *         asin_x = SUNRasin(x);
 * -----------------------------------------------------------------
 * SUNRasin(x) returns the asin value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRasin
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRasin(x) (asin((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRasin(x) (asinf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRasin(x) (asinl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRasin(x) (asinq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRacos
 * -----------------------------------------------------------------
 * Usage : sunrealtype acos_x;
 *         acos_x = SUNRacos(x);
 * -----------------------------------------------------------------
 * SUNRacos(x) returns the acos value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRacos
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRacos(x) (acos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRacos(x) (acosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRacos(x) (acosl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRacos(x) (acosq((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsinh
 * -----------------------------------------------------------------
 * Usage : sunrealtype sinh_x;
 *         sinh_x = SUNRsinh(x);
 * -----------------------------------------------------------------
 * SUNRsinh(x) returns sinh(x) (the hyperbolic sine of x).
 * -----------------------------------------------------------------
*/

#ifndef SUNRsinh
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsinh(x) (sinh((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsinh(x) (sinhf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsinh(x) (sinhl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRsinh(x) (sinhq((x)))
#else
#error \
"SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRcosh
 * -----------------------------------------------------------------
 * Usage : sunrealtype cosh_x;
 *         cosh_x = SUNRcosh(x);
 * -----------------------------------------------------------------
 * SUNRcosh(x) returns cosh(x) (the hyperbolic cosine of x).
 * -----------------------------------------------------------------
 */

#ifndef SUNRcosh
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRcosh(x) (cosh((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRcosh(x) (coshf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRcosh(x) (coshl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRcosh(x) (coshq((x)))
#else
#error \
"SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif


/*
 * -----------------------------------------------------------------
 * Function : SUNRatan
 * -----------------------------------------------------------------
 * Usage : sunrealtype cos_x;
 *         atan_x = SUNRatan(x);
 * -----------------------------------------------------------------
 * SUNRatan(x) returns the atan value of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNRatan
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRatan(x) (atan((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRatan(x) (atanf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRatan(x) (atanl((x)))
#elif defined(SUNDIALS_FLOAT128_PRECISION)
#define SUNRatan(x) (atanq((x)))
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
 * less than or equal to 10* or 100*(for __float128) machine epsilon.
 * If the relative difference is greater than 10* or 100*machine
 * epsilon, it returns 1. The function handles the case where a or b
 * are near zero as well as the case where a or b are inf/nan.
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
