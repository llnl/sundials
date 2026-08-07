
/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * Routines needed to read a sparse matrix in Rutherford-Boeing
 * format.
 * ----------------------------------------------------------------- */

#include <stdio.h>
#include <sundials/sundials_matrix.h>

#ifndef _DREADRB_H_
#define _DREADRB_H_

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

void dreadrb_dist(int iam, FILE* fp, SUNMatrix* Aout, SUNContext sunctx);

#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif

#endif