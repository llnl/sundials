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
 * Unit test for IDAS's internal sparse difference-quotient Jacobian routine.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "idas/idas.h"
#include "idas/idas_impl.h"
#include "idas/idas_ls.h"
#include "idas/idas_ls_impl.h"
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
  sunrealtype adata[NNZ];
  sunrealtype bdata[NNZ];
}* UserData;

static int check_retval(void* flagvalue, const char* funcname, int opt);
static int test_bad_sparse_template(SUNContext sunctx, UserData udata);
static int res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr,
               void* user_data);
static int fill_matrix(SUNMatrix A, UserData udata);
static sunrealtype matrix_value(sunindextype row, sunindextype col,
                                sunbooleantype yp_part);
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
  sunrealtype cj         = SUN_RCONST(2.0);
  SUNContext sunctx      = NULL;
  N_Vector y             = NULL;
  N_Vector yp            = NULL;
  N_Vector rr            = NULL;
  N_Vector tmp1          = NULL;
  N_Vector tmp2          = NULL;
  N_Vector tmp3          = NULL;
  SUNMatrix Jac          = NULL;
  SUNLinearSolver LS     = NULL;
  void* ida_mem          = NULL;
  IDAMem IDA_mem         = NULL;
  IDALsMem idals_mem     = NULL;
  sunindextype* groups   = NULL;
  sunindextype* rowmarks = NULL;

  struct UserData_ udata;

  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval(y, "N_VNew_Serial", 0)) return 1;
  yp = N_VClone(y);
  if (check_retval(yp, "N_VClone", 0)) return 1;
  rr = N_VClone(y);
  if (check_retval(rr, "N_VClone", 0)) return 1;
  tmp1 = N_VClone(y);
  if (check_retval(tmp1, "N_VClone", 0)) return 1;
  tmp2 = N_VClone(y);
  if (check_retval(tmp2, "N_VClone", 0)) return 1;
  tmp3 = N_VClone(y);
  if (check_retval(tmp3, "N_VClone", 0)) return 1;

  for (sunindextype i = 0; i < NEQ; i++)
  {
    NV_Ith_S(y, i)  = SUN_RCONST(1.0) + SUN_RCONST(0.1) * i;
    NV_Ith_S(yp, i) = -SUN_RCONST(0.5) + SUN_RCONST(0.05) * i;
  }

  Jac = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSC_MAT, sunctx);
  if (check_retval(Jac, "SUNSparseMatrix", 0)) return 1;

  retval = fill_matrix(Jac, &udata);
  if (check_retval(&retval, "fill_matrix", 1)) return 1;

  ida_mem = IDACreate(sunctx);
  if (check_retval(ida_mem, "IDACreate", 0)) return 1;

  retval = IDAInit(ida_mem, res, ZERO, y, yp);
  if (check_retval(&retval, "IDAInit", 1)) return 1;

  retval = IDASetUserData(ida_mem, &udata);
  if (check_retval(&retval, "IDASetUserData", 1)) return 1;

  retval = IDASStolerances(ida_mem, SUN_RCONST(1.0e-8), SUN_RCONST(1.0e-12));
  if (check_retval(&retval, "IDASStolerances", 1)) return 1;

  LS = test_linsol(sunctx);
  if (check_retval(LS, "test_linsol", 0)) return 1;

  retval = IDASetLinearSolver(ida_mem, LS, Jac);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) return 1;

  IDA_mem = (IDAMem)ida_mem;
  N_VConst(ONE, IDA_mem->ida_ewt);
  IDA_mem->ida_hh = ONE;

  retval = res(ZERO, y, yp, rr, &udata);
  if (check_retval(&retval, "res", 1)) return 1;

  retval = idaLsInitialize(IDA_mem);
  if (check_retval(&retval, "idaLsInitialize", 1)) return 1;

  idals_mem = (IDALsMem)IDA_mem->ida_lmem;
  if (check_retval(idals_mem, "ida_lmem", 0)) return 1;

  retval = idaLsDQJac(ZERO, cj, y, yp, rr, Jac, ida_mem, tmp1, tmp2, tmp3);
  if (check_retval(&retval, "idaLsDQJac", 1)) return 1;

  groups   = idals_mem->sparseDQgroups;
  rowmarks = idals_mem->sparseDQrowmarks;

  retval = idaLsDQJac(ZERO, cj, y, yp, rr, Jac, ida_mem, tmp1, tmp2, tmp3);
  if (check_retval(&retval, "idaLsDQJac", 1)) return 1;

  if ((groups != idals_mem->sparseDQgroups) ||
      (rowmarks != idals_mem->sparseDQrowmarks))
  {
    fprintf(stderr, "Sparse DQ Jacobian grouping workspace was not reused\n");
    fails++;
  }

  for (sunindextype j = 0; j < NEQ; j++)
  {
    for (sunindextype p = udata.colptrs[j]; p < udata.colptrs[j + 1]; p++)
    {
      sunrealtype exact = udata.adata[p] + cj * udata.bdata[p];
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

  if (idals_mem->nreDQ <= 0)
  {
    fprintf(stderr,
            "Expected internal sparse DQ Jacobian residual evaluations\n");
    fails++;
  }

  printf("Maximum sparse DQ Jacobian error = " SUN_FORMAT_G "\n", max_error);

  fails += test_bad_sparse_template(sunctx, &udata);

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(Jac);
  N_VDestroy(tmp3);
  N_VDestroy(tmp2);
  N_VDestroy(tmp1);
  N_VDestroy(rr);
  N_VDestroy(yp);
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
  N_Vector yp        = NULL;
  SUNMatrix Jac      = NULL;
  SUNLinearSolver LS = NULL;
  SUNLogger logger   = NULL;
  FILE* errfp        = NULL;
  void* ida_mem      = NULL;
  IDAMem IDA_mem     = NULL;

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval(y, "N_VNew_Serial", 0)) return 1;

  yp = N_VClone(y);
  if (check_retval(yp, "N_VClone", 0)) return 1;

  Jac = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSR_MAT, sunctx);
  if (check_retval(Jac, "SUNSparseMatrix", 0)) return 1;

  LS = test_linsol(sunctx);
  if (check_retval(LS, "test_linsol", 0)) return 1;

  ida_mem = IDACreate(sunctx);
  if (check_retval(ida_mem, "IDACreate", 0)) return 1;

  retval = IDAInit(ida_mem, res, ZERO, y, yp);
  if (check_retval(&retval, "IDAInit", 1)) fails++;

  retval = IDASetUserData(ida_mem, udata);
  if (check_retval(&retval, "IDASetUserData", 1)) fails++;

  retval = IDASetLinearSolver(ida_mem, LS, Jac);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) fails++;

  IDA_mem = (IDAMem)ida_mem;

  retval = SUNContext_GetLogger(sunctx, &logger);
  if (check_retval(&retval, "SUNContext_GetLogger", 1)) fails++;

  retval = SUNLogger_GetErrorFile(logger, &errfp);
  if (check_retval(&retval, "SUNLogger_GetErrorFile", 1)) fails++;

  retval = SUNLogger_SetErrorFile(logger, NULL);
  if (check_retval(&retval, "SUNLogger_SetErrorFile", 1)) fails++;

  init_retval = idaLsInitialize(IDA_mem);

  retval = SUNLogger_SetErrorFile(logger, errfp);
  if (check_retval(&retval, "SUNLogger_SetErrorFile", 1)) fails++;

  if (init_retval != IDALS_ILL_INPUT)
  {
    fprintf(stderr,
            "Expected idaLsInitialize to reject CSR sparse DQ template with "
            "IDALS_ILL_INPUT, retval = %d\n",
            init_retval);
    fails++;
  }

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(Jac);
  N_VDestroy(yp);
  N_VDestroy(y);

  return fails;
}

static int res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr,
               void* user_data)
{
  UserData udata      = (UserData)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* ypdata = N_VGetArrayPointer(yp);
  sunrealtype* rdata  = N_VGetArrayPointer(rr);
  sunindextype j      = 0;
  sunindextype p      = 0;

  N_VConst(ZERO, rr);

  for (j = 0; j < NEQ; j++)
  {
    for (p = udata->colptrs[j]; p < udata->colptrs[j + 1]; p++)
    {
      rdata[udata->rowvals[p]] += udata->adata[p] * ydata[j] +
                                  udata->bdata[p] * ypdata[j];
    }
  }

  return 0;
}

static int fill_matrix(SUNMatrix A, UserData udata)
{
  sunindextype p = 0;

  if (SUNSparseMatrix_SparseType(A) != SUN_CSC_MAT) { return -1; }

  for (sunindextype j = 0; j < NEQ; j++)
  {
    udata->colptrs[j] = p;

    udata->rowvals[p] = j;
    udata->adata[p]   = matrix_value(j, j, SUNFALSE);
    udata->bdata[p]   = matrix_value(j, j, SUNTRUE);
    p++;

    udata->rowvals[p] = (j + 3) % NEQ;
    udata->adata[p]   = matrix_value((j + 3) % NEQ, j, SUNFALSE);
    udata->bdata[p]   = matrix_value((j + 3) % NEQ, j, SUNTRUE);
    p++;

    udata->rowvals[p] = (j + 7) % NEQ;
    udata->adata[p]   = matrix_value((j + 7) % NEQ, j, SUNFALSE);
    udata->bdata[p]   = matrix_value((j + 7) % NEQ, j, SUNTRUE);
    p++;
  }
  udata->colptrs[NEQ] = p;

  for (p = 0; p < NEQ + 1; p++) { SM_INDEXPTRS_S(A)[p] = udata->colptrs[p]; }
  for (p = 0; p < NNZ; p++)
  {
    SM_INDEXVALS_S(A)[p] = udata->rowvals[p];
    SM_DATA_S(A)[p]      = ZERO;
  }

  return (p == NNZ) ? 0 : -1;
}

static sunrealtype matrix_value(sunindextype row, sunindextype col,
                                sunbooleantype yp_part)
{
  sunrealtype scale = yp_part ? SUN_RCONST(0.0625) : SUN_RCONST(0.125);

  if (row == col)
  {
    return (yp_part ? SUN_RCONST(0.5) : -SUN_RCONST(2.0)) - scale * (col + 1);
  }
  return scale * (row + 1) + SUN_RCONST(0.25) * (col + 1);
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
