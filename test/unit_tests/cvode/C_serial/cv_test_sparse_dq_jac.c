/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * Unit test for CVODE's internal sparse difference-quotient Jacobian routine.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "cvode/cvode.h"
#include "cvode/cvode_impl.h"
#include "cvode/cvode_ls.h"
#include "cvode/cvode_ls_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_linearsolver.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_math.h"
#include "sunmatrix/sunmatrix_sparse.h"

#define NEQ 10
#define NNZ 30

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

typedef struct UserData_
{
  sunindextype colptrs[NEQ + 1];
  sunindextype rowvals[NNZ];
  sunrealtype data[NNZ];
}* UserData;

static int check_retval(void* flagvalue, const char* funcname, int opt);
static int test_bad_sparse_template(SUNContext sunctx, UserData udata);
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fill_matrix(SUNMatrix A, UserData udata, sunbooleantype copy_data);
static sunrealtype matrix_value(sunindextype row, sunindextype col);
static SUNLinearSolver test_linsol(SUNContext sunctx);
static SUNLinearSolver_Type test_linsol_gettype(SUNLinearSolver S);
static int test_linsol_solve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, sunrealtype tol);
static SUNErrCode test_linsol_free(SUNLinearSolver S);

int main(void)
{
  int fails              = 0;
  int retval             = 0;
  sunrealtype max_error  = ZERO;
  SUNContext sunctx      = NULL;
  N_Vector y             = NULL;
  N_Vector fy            = NULL;
  N_Vector tmp1          = NULL;
  N_Vector tmp2          = NULL;
  N_Vector tmp3          = NULL;
  SUNMatrix Jac          = NULL;
  SUNLinearSolver LS     = NULL;
  void* cvode_mem        = NULL;
  CVodeMem cv_mem        = NULL;
  CVLsMem cvls_mem       = NULL;
  sunindextype* groups   = NULL;
  sunindextype* rowmarks = NULL;

  struct UserData_ udata;

  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval(y, "N_VNew_Serial", 0)) return 1;
  fy = N_VClone(y);
  if (check_retval(fy, "N_VClone", 0)) return 1;
  tmp1 = N_VClone(y);
  if (check_retval(tmp1, "N_VClone", 0)) return 1;
  tmp2 = N_VClone(y);
  if (check_retval(tmp2, "N_VClone", 0)) return 1;
  tmp3 = N_VClone(y);
  if (check_retval(tmp3, "N_VClone", 0)) return 1;

  Jac = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSC_MAT, sunctx);
  if (check_retval(Jac, "SUNSparseMatrix", 0)) return 1;

  retval = fill_matrix(Jac, &udata, SUNFALSE);
  if (check_retval(&retval, "fill_matrix", 1)) return 1;

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval(cvode_mem, "CVodeCreate", 0)) return 1;

  retval = CVodeInit(cvode_mem, rhs, ZERO, y);
  if (check_retval(&retval, "CVodeInit", 1)) return 1;

  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return 1;

  retval = CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-8), SUN_RCONST(1.0e-12));
  if (check_retval(&retval, "CVodeSStolerances", 1)) return 1;

  LS = test_linsol(sunctx);
  if (check_retval(LS, "test_linsol", 0)) return 1;

  retval = CVodeSetLinearSolver(cvode_mem, LS, Jac);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

  cv_mem = (CVodeMem)cvode_mem;
  N_VConst(ONE, cv_mem->cv_ewt);
  cv_mem->cv_h = ONE;

  retval = rhs(ZERO, y, fy, &udata);
  if (check_retval(&retval, "rhs", 1)) return 1;

  retval = cvLsInitialize(cv_mem);
  if (check_retval(&retval, "cvLsInitialize", 1)) return 1;

  cvls_mem = (CVLsMem)cv_mem->cv_lmem;
  if (check_retval(cvls_mem, "cv_lmem", 0)) return 1;

  retval = cvLsDQJac(ZERO, y, fy, Jac, cvode_mem, tmp1, tmp2, tmp3);
  if (check_retval(&retval, "cvLsDQJac", 1)) return 1;

  groups   = cvls_mem->sparseDQgroups;
  rowmarks = cvls_mem->sparseDQrowmarks;

  retval = cvLsDQJac(ZERO, y, fy, Jac, cvode_mem, tmp1, tmp2, tmp3);
  if (check_retval(&retval, "cvLsDQJac", 1)) return 1;

  if ((groups != cvls_mem->sparseDQgroups) ||
      (rowmarks != cvls_mem->sparseDQrowmarks))
  {
    fprintf(stderr, "Sparse DQ Jacobian grouping workspace was not reused\n");
    fails++;
  }

  for (sunindextype j = 0; j < NEQ; j++)
  {
    for (sunindextype p = udata.colptrs[j]; p < udata.colptrs[j + 1]; p++)
    {
      sunrealtype exact = udata.data[p];
      sunrealtype error = SUNRabs(SM_DATA_S(Jac)[p] - exact);
      max_error         = SUNMAX(max_error, error);
      if (error > SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF)
      {
        fprintf(stderr,
                "Jacobian mismatch at row %lld, column %lld: computed "
                "= " SUN_FORMAT_G ", expected = " SUN_FORMAT_G
                ", error = " SUN_FORMAT_G "\n",
                (long long)udata.rowvals[p], (long long)j, SM_DATA_S(Jac)[p],
                exact, error);
        fails++;
      }
    }
  }

  if (cvls_mem->nfeDQ <= 0)
  {
    fprintf(stderr, "Expected internal sparse DQ Jacobian RHS evaluations\n");
    fails++;
  }

  printf("Maximum sparse DQ Jacobian error = " SUN_FORMAT_G "\n", max_error);

  fails += test_bad_sparse_template(sunctx, &udata);

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(Jac);
  N_VDestroy(tmp3);
  N_VDestroy(tmp2);
  N_VDestroy(tmp1);
  N_VDestroy(fy);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return fails;
}

static int test_bad_sparse_template(SUNContext sunctx, UserData udata)
{
  int fails          = 0;
  int init_retval    = 0;
  int retval         = 0;
  N_Vector y         = NULL;
  SUNMatrix Jac      = NULL;
  SUNLinearSolver LS = NULL;
  SUNLogger logger   = NULL;
  FILE* errfp        = NULL;
  void* cvode_mem    = NULL;
  CVodeMem cv_mem    = NULL;

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval(y, "N_VNew_Serial", 0)) return 1;

  Jac = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSR_MAT, sunctx);
  if (check_retval(Jac, "SUNSparseMatrix", 0)) return 1;

  LS = test_linsol(sunctx);
  if (check_retval(LS, "test_linsol", 0)) return 1;

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval(cvode_mem, "CVodeCreate", 0)) return 1;

  retval = CVodeInit(cvode_mem, rhs, ZERO, y);
  if (check_retval(&retval, "CVodeInit", 1)) fails++;

  retval = CVodeSetUserData(cvode_mem, udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) fails++;

  retval = CVodeSetLinearSolver(cvode_mem, LS, Jac);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) fails++;

  cv_mem = (CVodeMem)cvode_mem;

  retval = SUNContext_GetLogger(sunctx, &logger);
  if (check_retval(&retval, "SUNContext_GetLogger", 1)) fails++;

  retval = SUNLogger_GetErrorFile(logger, &errfp);
  if (check_retval(&retval, "SUNLogger_GetErrorFile", 1)) fails++;

  retval = SUNLogger_SetErrorFile(logger, NULL);
  if (check_retval(&retval, "SUNLogger_SetErrorFile", 1)) fails++;

  init_retval = cvLsInitialize(cv_mem);

  retval = SUNLogger_SetErrorFile(logger, errfp);
  if (check_retval(&retval, "SUNLogger_SetErrorFile", 1)) fails++;

  if (init_retval != CVLS_ILL_INPUT)
  {
    fprintf(stderr,
            "Expected cvLsInitialize to reject CSR sparse DQ template with "
            "CVLS_ILL_INPUT, retval = %d\n",
            init_retval);
    fails++;
  }

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(Jac);
  N_VDestroy(y);

  return fails;
}

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData udata     = (UserData)user_data;
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);
  sunindextype j     = 0;
  sunindextype p     = 0;

  N_VConst(ZERO, ydot);

  for (j = 0; j < NEQ; j++)
  {
    for (p = udata->colptrs[j]; p < udata->colptrs[j + 1]; p++)
    {
      fdata[udata->rowvals[p]] += udata->data[p] * ydata[j];
    }
  }

  return 0;
}

static int fill_matrix(SUNMatrix A, UserData udata, sunbooleantype copy_data)
{
  sunindextype p = 0;

  if (SUNSparseMatrix_SparseType(A) != SUN_CSC_MAT) { return -1; }

  for (sunindextype j = 0; j < NEQ; j++)
  {
    udata->colptrs[j] = p;

    udata->rowvals[p] = j;
    udata->data[p]    = matrix_value(j, j);
    p++;

    udata->rowvals[p] = (j + 3) % NEQ;
    udata->data[p]    = matrix_value((j + 3) % NEQ, j);
    p++;

    udata->rowvals[p] = (j + 7) % NEQ;
    udata->data[p]    = matrix_value((j + 7) % NEQ, j);
    p++;
  }
  udata->colptrs[NEQ] = p;

  for (p = 0; p < NEQ + 1; p++) { SM_INDEXPTRS_S(A)[p] = udata->colptrs[p]; }
  for (p = 0; p < NNZ; p++)
  {
    SM_INDEXVALS_S(A)[p] = udata->rowvals[p];
    SM_DATA_S(A)[p]      = copy_data ? udata->data[p] : ZERO;
  }

  return (p == NNZ) ? 0 : -1;
}

static sunrealtype matrix_value(sunindextype row, sunindextype col)
{
  if (row == col) { return -SUN_RCONST(2.0) - SUN_RCONST(0.25) * (col + 1); }
  return SUN_RCONST(0.125) * (row + 1) + SUN_RCONST(0.25) * (col + 1);
}

static SUNLinearSolver test_linsol(SUNContext sunctx)
{
  SUNLinearSolver LS = SUNLinSolNewEmpty(sunctx);
  if (LS == NULL) { return NULL; }

  LS->ops->gettype = test_linsol_gettype;
  LS->ops->solve   = test_linsol_solve;
  LS->ops->free    = test_linsol_free;

  return LS;
}

static SUNLinearSolver_Type test_linsol_gettype(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

static int test_linsol_solve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, sunrealtype tol)
{
  N_VScale(ONE, b, x);
  return SUN_SUCCESS;
}

static SUNErrCode test_linsol_free(SUNLinearSolver S)
{
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a value so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
