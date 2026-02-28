/*
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * This is the header file contains the prototypes for functions to
 * test SUNMatrix module implementation.
 * -----------------------------------------------------------------
 */

#include <math.h>
#include <sundials/sundials_core.h>

/* define constants */
#define NEG_TWO  SUN_RCONST(-2.0)
#define NEG_ONE  SUN_RCONST(-1.0)
#define NEG_HALF SUN_RCONST(-0.5)
#define ZERO     SUN_RCONST(0.0)
#define HALF     SUN_RCONST(0.5)
#define ONE      SUN_RCONST(1.0)
#define TWO      SUN_RCONST(2.0)
#define THREE    SUN_RCONST(3.0)

/* Helpers for printing out test status information */
#define TEST_STATUS(fmt, myrank)                            \
  if (print_all_ranks == 0 && myrank == 0) { printf(fmt); } \
  else if (print_all_ranks)                                 \
  {                                                         \
    printf("process %06d: ", myrank);                       \
    printf(fmt);                                            \
  }
#define TEST_STATUS2(fmt, retval, myrank)                           \
  if (print_all_ranks == 0 && myrank == 0) { printf(fmt, retval); } \
  else if (print_all_ranks)                                         \
  {                                                                 \
    printf("process %06d: ", myrank);                               \
    printf(fmt, retval);                                            \
  }
#define TEST_STATUS3(fmt, val1, val2, myrank)                           \
  if (print_all_ranks == 0 && myrank == 0) { printf(fmt, val1, val2); } \
  else if (print_all_ranks)                                             \
  {                                                                     \
    printf("process %06d: ", myrank);                                   \
    printf(fmt, val1, val2);                                            \
  }

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif
/* Forward declarations for implementation specific utility functions */
int check_matrix(SUNMatrix A, SUNMatrix B, sunrealtype tol);
int check_matrix_entry(SUNMatrix A, sunrealtype val, sunrealtype tol);
int check_vector(N_Vector expected, N_Vector computed, sunrealtype tol);
sunbooleantype has_data(SUNMatrix A);
sunbooleantype is_square(SUNMatrix A);
void sync_device(SUNMatrix A);

/* Test function declarations */
int Test_SUNMatGetID(SUNMatrix A, SUNMatrix_ID sunid, int myid);
int Test_SUNMatClone(SUNMatrix A, int myid);
int Test_SUNMatZero(SUNMatrix A, int myid);
int Test_SUNMatCopy(SUNMatrix A, int myid);
int Test_SUNMatScaleAdd(SUNMatrix A, SUNMatrix I, int myid);
int Test_SUNMatScaleAddI(SUNMatrix A, SUNMatrix I, int myid);
int Test_SUNMatMatvecSetup(SUNMatrix A, int myid);
int Test_SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y, int myid);
int Test_SUNMatHermitianTransposeVec(SUNMatrix A, SUNMatrix AT, N_Vector x,
                                     N_Vector y, int myid);
int Test_SUNMatSpace(SUNMatrix A, int myid);

/* Timing function */
void SetTiming(int onoff);

/* Function to set print_all_ranks */
void SetPrintAllRanks(int onoff);
#ifdef __cplusplus
}
#endif
