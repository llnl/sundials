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
 * Macro : SUN_COMPLEX
 * -----------------------------------------------------------------
 * Usage : suncomplextype x = SUN_COMPLEX(ONE, ONE);
 * -----------------------------------------------------------------
 * SUN_COMPLEX(real, imag) returns a complex number with real part
 * real and imaginary part imag. This macro is used to initialize
 * complex numbers.
 * -----------------------------------------------------------------
 */

#ifndef SUN_COMPLEX
#define SUN_COMPLEX(real, imag) (real + imag * SUN_I)
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUN_CREAL, SUN_CIMAG, SUN_REAL, SUN_IMAG
 * -----------------------------------------------------------------
 * Usage : sunscalartype x;
 *         suncomplextype z;
 *         sunrealtype real_x, imag_x, real_z, imag_z;
 *         real_z = SUN_CREAL(z);
 *         imag_z = SUN_CIMAG(z);
 *         real_x = SUN_REAL(x);
 *         real_x = SUN_REAL(x);
 * -----------------------------------------------------------------
 * SUNCREAL and SUNCIMAG return the real and imaginary components of
 *   their suncomplextype arguments, respectively.
 * If sunscalartype is complex, then SUNREAL and SUNIMAG just call
 *   SUNCREAL and SUNCIMAG.
 * If sunscalartype is real, then SUNREAL just returns its argument
 *   and SUNIMAG returns zero.
 * -----------------------------------------------------------------
 */

#ifndef SUN_CREAL
#ifdef __cplusplus
#define SUN_CREAL(z) (std::real(z))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUN_CREAL(z) (creal(z))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUN_CREAL(z) (crealf(z))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUN_CREAL(z) (creall(z))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

#ifndef SUN_CIMAG
#ifdef __cplusplus
#define SUN_CIMAG(z) (std::imag(z))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUN_CIMAG(z) (cimag(z))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUN_CIMAG(z) (cimagf(z))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUN_CIMAG(z) (cimagl(z))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUN_REAL(z) (z)
#define SUN_IMAG(z) SUN_RCONST(0.0)
#elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
#define SUN_REAL(z) SUN_CREAL(z)
#define SUN_IMAG(z) SUN_CIMAG(z)
#else
#error \
  "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNCCONJ
 * -----------------------------------------------------------------
 * Usage : suncomplextype sqrt_x;
 *         sqrt_x = SUNCCONJ(x);
 * -----------------------------------------------------------------
 * SUNCCONJ(x) returns the complex conjugate of x.
 * -----------------------------------------------------------------
 */

#ifndef SUNCCONJ
#ifdef __cplusplus
#define SUNCCONJ(z) (std::conj(z))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCCONJ(z) (conj(z))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCCONJ(z) (conjf(z))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCCONJ(z) (conjl(z))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNCONJ
 * -----------------------------------------------------------------
 * Usage : sunscalartype sqrt_x;
 *         sqrt_x = SUNCONJ(x);
 * -----------------------------------------------------------------
 * SUNCONJ(x) returns the complex conjugate of x if x is complex;
 * if x is real then it does nothing.
 * -----------------------------------------------------------------
 */

#ifndef SUNCONJ
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNCONJ(x) (x)
#else
#define SUNCONJ(x) SUNCCONJ(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Macros
 * -----------------------------------------------------------------
 * SUNMIN(A,B) returns the minimum of A and B
 * SUNCMIN(A,B) returns whichever of A and B has minimum real part
 * SUNSMIN(A,B) uses SUNMIN if sunscalartype is real, and SUNCMIN if
 *   sunscalartype is complex
 *
 * SUNMAX(A,B) returns the maximum of A and B
 * SUNCMAX(A,B) returns whichever of A and B has maximum real part
 * SUNSMAX(A,B) uses SUNMAX if sunscalartype is real, and SUNCMAX if
 *   sunscalartype is complex
 *
 * SUNRSQR(A) returns A^2
 * SUNCSQR(A) returns A*conj(A)
 * SUNSQR(A) uses whichever of the above is mapped to sunscalartype
 *
 * SUNRsqrt calls the appropriate version of sqrt (real numbers)
 * SUNCsqrt calls the appropriate version of csqrt
 * SUNsqrt uses whichever of the above is mapped to sunrealtype
 *
 * SUNRabs calls the appropriate version of abs (real numbers)
 * SUNCabs calls the appropriate version of cabs
 * SUNabs uses whichever of the above is mapped to sunrealtype
 *
 * SUNRexp calls the appropriate version of exp (real numbers)
 * SUNCexp calls the appropriate version of cexp
 * SUNexp uses whichever of the above is mapped to sunrealtype
 *
 * SUNRceil calls the appropriate version of ceil (real numbers)
 *
 * SUNRround calls the appropriate version of round
 * -----------------------------------------------------------------
 */

#ifndef SUNMIN
#define SUNMIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef SUNCMIN
#define SUNCMIN(A, B) ((SUN_REAL(A)) < (SUN_REAL(B)) ? (A) : (B))
#endif

#ifndef SUNSMIN
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNSMIN(A, B) SUNMIN(A, B)
#else
#define SUNSMIN(A, B) SUNCMIN(A, B)
#endif
#endif

#ifndef SUNMAX
#define SUNMAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef SUNCMAX
#define SUNCMAX(A, B) ((SUN_REAL(A)) > (SUN_REAL(B)) ? (A) : (B))
#endif

#ifndef SUNSMAX
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNSMAX(A, B) SUNMAX(A, B)
#else
#define SUNSMAX(A, B) SUNCMAX(A, B)
#endif
#endif

#ifndef SUNRSQR
#define SUNRSQR(A) ((A) * (A))
#endif

#ifndef SUNCSQR
#define SUNCSQR(A) SUN_REAL((A) * SUNCCONJ(A))
#endif

#ifndef SUNSQR
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNSQR(A) SUNRSQR(A)
#else
#define SUNSQR(A) SUN_REAL(SUNCSQR(A))
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsqrt, SUNCsqrt, and SUNsqrt
 * -----------------------------------------------------------------
 * Usage : sunrealtype sqrt_x;
 *         sqrt_x = SUNRsqrt(x);
 *         suncomplextype sqrt_y;
 *         sqrt_y = SUNCsqrt(y);
 *         sunscalartype sqrt_z;
 *         sqrt_z = SUNsqrt(z);
 * -----------------------------------------------------------------
 * SUNRsqrt(x) returns the square root of x. If x < ZERO, then
 * SUNRsqrt returns ZERO.
 * SUNCsqrt(y) returns the complex square root of x.
 * SUNsqrt(z) uses whichever of the above is mapped to sunscalartype
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

#ifndef SUNCsqrt
#ifdef __cplusplus
#define SUNCsqrt(z) (std::sqrt(z))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCsqrt(z) (csqrt(z))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCsqrt(z) (csqrtf(z))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCsqrt(z) (csqrtl(z))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

#ifndef SUNsqrt
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNsqrt(x) SUNRsqrt(x)
#else
#define SUNsqrt(x) SUNCsqrt(x)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRabs, SUNCabs, SUNabs
 * -----------------------------------------------------------------
 * Usage : sunrealtype abs_x;
 *         abs_x = SUNRabs(x);
 *         suncomplextype abs_y;
 *         abs_y = SUNCabs(y);
 *         sunscalartype abs_z;
 *         abs_z = SUNabs(z);
 * -----------------------------------------------------------------
 * SUNRabs(x) returns the absolute value of x.
 * SUNCabs(x) returns the complex absolute value of x.
 * SUNabs(x) uses whichever of the above is mapped to sunscalartype.
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

#ifndef SUNCabs
#ifdef __cplusplus
#define SUNCabs(z) (std::abs(z))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCabs(z) (cabs(z))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCabs(z) (cabsf(z))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCabs(z) (cabsl(z))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

#ifndef SUNabs
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNabs(x) SUNRabs(x)
#else
#define SUNabs(x) SUN_REAL(SUNCabs(x))
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRexp, SUNCexp, SUNexp
 * -----------------------------------------------------------------
 * Usage : sunrealtype exp_x;
 *         exp_x = SUNRexp(x);
 *         suncomplextype exp_y;
 *         exp_y = SUNCexp(y);
 *         sunscalartype exp_z;
 *         exp_z = SUNexp(z);
 * -----------------------------------------------------------------
 * SUNRexp(x) returns e^x (base-e exponential function) for
 * real-valued x.
 * SUNCexp(y) returns e^y (base-e exponential function) for
 * complex-valued y.
 * SUNexp(z) uses whichever of the above is mapped to sunscalartype.
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

#ifndef SUNCexp
#ifdef __cplusplus
#define SUNCexp(z) (std::exp(z))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCexp(z) (cexp(z))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCexp(z) (cexpf(z))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCexp(z) (cexpl(z))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

#ifndef SUNexp
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNexp(x) SUNRexp(x)
#else
#define SUNexp(x) SUNCexp(x)
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
 *         sunrealtype rbase, rans;
 *         rans = SUNRpowerI(rbase,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerI returns the value of base^exponent, where base is of
 * type sunrealtype and exponent is of type int.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunrealtype SUNRpowerI(sunrealtype base, int exponent);

/*
 * -----------------------------------------------------------------
 * Function : SUNRpowerR, SUNCpowerC, SUNpower
 * -----------------------------------------------------------------
 * Usage : sunrealtype rbase, rexponent, rans;
 *         rans = SUNRpowerR(rbase,rexponent);
 *         suncomplextype cbase, cexponent, cans;
 *         cans = SUNRpowerR(cbase,cexponent);
 *         sunscalartype base, exponent, ans;
 *         ans = SUNpower(base,exponent);
 * -----------------------------------------------------------------
 * SUNRpowerR returns the value of base^exponent, where both base and
 * exponent are of type sunrealtype. If base < ZERO, then SUNRpowerR
 * returns ZERO.
 * SUNCpowerC returns the value of base^exponent, where both base and
 * exponent are of type suncomplextype.
 * SUNpower uses whichever of the above is mapped to sunscalartype.
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

#ifndef SUNCpowerC
#ifdef __cplusplus
#define SUNCpowerC(base, exponent) (std::pow((base), (exponent)))
#else
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNCpowerC(base, exponent) ((cpow((base), (exponent))))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNCpowerC(base, exponent) ((cpowf((base), (exponent))))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNCpowerC(base, exponent) ((cpowl((base), (exponent))))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif
#endif

#ifndef SUNpower
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNpower(base, exponent) SUNRpowerR(base, exponent)
#else
#define SUNpower(base, exponent) SUNCpowerC(base, exponent)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompare, SUNCCompare, SUNCompare
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype ra, rb;
 *         isNotEqual = SUNRCompare(ra, rb);
 *         suncomplextype ca, cb;
 *         isNotEqual = SUNCCompare(ca, cb);
 *         sunscalartype a, b;
 *         isNotEqual = SUNCompare(a, b);
 * -----------------------------------------------------------------
 * These functions return 0 if the relative difference of a and b is
 * less than or equal to 10*machine epsilon. If the relative
 * difference is greater than 10*machine epsilon, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * SUNRCompare is designed for real-valued inputs.
 * SUNCCompare is designed for complex-valued inputs.
 * SUNCompare uses whichever of the above is mapped to sunscalartype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);
SUNDIALS_EXPORT sunbooleantype SUNCCompare(suncomplextype a, suncomplextype b);

#ifndef SUNCompare
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNCompare(a, b) SUNRCompare(a, b)
#else
#define SUNCompare(a, b) SUNCCompare(a, b)
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRCompareTol, SUNCompareTol, SUNCompareTol
 * -----------------------------------------------------------------
 * Usage : int isNotEqual;
 *         sunrealtype tol;
 *         sunrealtype ra, rb;
 *         isNotEqual = SUNRCompareTol(ra, rb, tol);
 *         suncomplextype ca, cb;
 *         isNotEqual = SUNRCompareTol(ca, cb, tol);
 *         sunscalartype a, b;
 *         isNotEqual = SUNCompareTol(a, b, tol);
 * -----------------------------------------------------------------
 * These return 0 if the relative difference of a and b is
 * less than or equal to the provided tolerance. If the relative
 * difference is greater than the tolerance, it returns 1. The
 * function handles the case where a or b are near zero as well as
 * the case where a or b are inf/nan.
 * SUNRCompareTol is designed for real-valued inputs.
 * SUNCCompareTol is designed for complex-valued inputs.
 * SUNCompareTol uses whichever of the above is mapped to
 * sunscalartype.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b,
                                              sunrealtype tol);
SUNDIALS_EXPORT sunbooleantype SUNCCompareTol(suncomplextype a,
                                              suncomplextype b, sunrealtype tol);

#ifndef SUNCompareTol
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define SUNCompareTol(a, b, tol) SUNRCompareTol(a, b, tol)
#else
#define SUNCompareTol(a, b, tol) SUNCCompareTol(a, b, tol)
#endif
#endif

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
