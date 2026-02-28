/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL, Daniel R. Reynolds @ UMBC
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
 * These are test functions for an NVECTOR module implementation
 * which have MPI symbols.
 * -----------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_nvector.h>

#include "test_nvector.h"

void Test_AbortMPI(SUNComm comm, int code)
{
  Test_Finalize();
  MPI_Abort(comm, code);
}

/* ----------------------------------------------------------------------
 * Test_N_VGetCommunicator Test (with MPI dependency).
 * --------------------------------------------------------------------*/
int Test_N_VGetCommunicatorMPI(N_Vector W, SUNComm comm, int myid)
{
  SUNComm wcomm;
  int same;

  /* ask W for its communicator */
  wcomm = N_VGetCommunicator(W);

  /* return with success if both are NULL */
  if ((wcomm == SUN_COMM_NULL) && (comm == SUN_COMM_NULL))
  {
    printf("PASSED test -- N_VGetCommunicator\n");
    return (0);
  }

  /* return with failure if either is NULL */
  if (wcomm == SUN_COMM_NULL)
  {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (incorrectly "
           "reports NULL comm)\n",
           myid);
    return (1);
  }
  if (comm == SUN_COMM_NULL)
  {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (incorrectly "
           "reports non-NULL comm)\n",
           myid);
    return (1);
  }

  /* call MPI_Comm_compare to check that communicators match or are congruent */
  if (MPI_Comm_compare(comm, wcomm, &same) != MPI_SUCCESS)
  {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (error in "
           "MPI_Comm_compare)\n",
           myid);
    return (1);
  }
  if ((same != MPI_IDENT) && (same != MPI_CONGRUENT))
  {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (mismatched "
           "comms)\n",
           myid);
    return (1);
  }
  if (myid == 0) { printf("PASSED test -- N_VGetCommunicator\n"); }
  return (0);
}
