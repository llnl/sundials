/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ UMBC
 * ---------------------------------------------------------------------------
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
 * ---------------------------------------------------------------------------
 * This is the implementation file for ARKODE's MRIStepCoupling tables.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"

/* ===========================================================================
 * Exported Functions
 * ===========================================================================*/

/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  imeth -- integer key for the desired method
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID method)
{
  switch (method)
  {
#define ARK_MRI_TABLE(name, coeff) \
  case name: coeff break;
#include "arkode_mri_tables.def"
#undef ARK_MRI_TABLE

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown coupling table");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  method -- string key for the desired method
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_LoadTableByName(const char* method)
{
#define ARK_MRI_TABLE(name, coeff) \
  if (strcmp(#name, method) == 0) coeff
#include "arkode_mri_tables.def"
#undef ARK_MRI_TABLE

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown coupling table");

  return NULL;
}

/*---------------------------------------------------------------
  Routine to allocate an empty MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages,
                                      MRISTEP_METHOD_TYPE type)
{
  int i, j;
  sunbooleantype hasOmegas, hasGammas;
  MRIStepCoupling MRIC = NULL;

  /* Check for legal input values */
  if (nmat < 1 || stages < 1) { return (NULL); }

  /* ------------------------------------------
   * Allocate and initialize coupling structure
   * ------------------------------------------ */

  MRIC = (MRIStepCoupling)malloc(sizeof(struct MRIStepCouplingMem));
  if (!MRIC) { return (NULL); }

  MRIC->type   = type;
  MRIC->nmat   = nmat;
  MRIC->stages = stages;
  MRIC->q      = 0;
  MRIC->p      = 0;
  MRIC->c      = NULL;
  MRIC->W      = NULL;
  MRIC->G      = NULL;
  MRIC->ngroup = 0;
  MRIC->group  = NULL;

  /* --------------------------------------------
   * Determine general storage format
   * -------------------------------------------- */

  hasOmegas = hasGammas = SUNFALSE;
  if ((type == MRISTEP_EXPLICIT) || (type == MRISTEP_IMEX) ||
      (type == MRISTEP_MERK) || (type == MRISTEP_SR))
  {
    hasOmegas = SUNTRUE;
  }
  if ((type == MRISTEP_IMPLICIT) || (type == MRISTEP_IMEX) ||
      (type == MRISTEP_SR))
  {
    hasGammas = SUNTRUE;
  }

  /* --------------------------------------------
   * Allocate abscissae and coupling coefficients
   * -------------------------------------------- */

  MRIC->c = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
  if (!(MRIC->c))
  {
    MRIStepCoupling_Free(MRIC);
    return (NULL);
  }

  if (hasOmegas)
  {
    /* allocate W matrices */
    MRIC->W = (sunrealtype***)calloc(nmat, sizeof(sunrealtype**));
    if (!(MRIC->W))
    {
      MRIStepCoupling_Free(MRIC);
      return (NULL);
    }

    /* allocate rows of each matrix in W */
    for (i = 0; i < nmat; i++)
    {
      MRIC->W[i] = NULL;
      MRIC->W[i] = (sunrealtype**)calloc(stages + 1, sizeof(sunrealtype*));
      if (!(MRIC->W[i]))
      {
        MRIStepCoupling_Free(MRIC);
        return (NULL);
      }
    }

    /* allocate columns of each matrix in W */
    for (i = 0; i < nmat; i++)
    {
      for (j = 0; j <= stages; j++)
      {
        MRIC->W[i][j] = NULL;
        MRIC->W[i][j] = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
        if (!(MRIC->W[i][j]))
        {
          MRIStepCoupling_Free(MRIC);
          return (NULL);
        }
      }
    }
  }

  if (hasGammas)
  {
    /* allocate G matrices */
    MRIC->G = (sunrealtype***)calloc(nmat, sizeof(sunrealtype**));
    if (!(MRIC->G))
    {
      MRIStepCoupling_Free(MRIC);
      return (NULL);
    }

    /* allocate rows of each matrix in G */
    for (i = 0; i < nmat; i++)
    {
      MRIC->G[i] = NULL;
      MRIC->G[i] = (sunrealtype**)calloc(stages + 1, sizeof(sunrealtype*));
      if (!(MRIC->G[i]))
      {
        MRIStepCoupling_Free(MRIC);
        return (NULL);
      }
    }

    /* allocate columns of each matrix in G */
    for (i = 0; i < nmat; i++)
    {
      for (j = 0; j <= stages; j++)
      {
        MRIC->G[i][j] = NULL;
        MRIC->G[i][j] = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
        if (!(MRIC->G[i][j]))
        {
          MRIStepCoupling_Free(MRIC);
          return (NULL);
        }
      }
    }
  }

  /* for MERK methods, allocate maximum possible number/sizes of stage groups */
  if (type == MRISTEP_MERK)
  {
    MRIC->ngroup = stages;
    MRIC->group  = (int**)malloc(stages * sizeof(int*));
    if (!(MRIC->group))
    {
      MRIStepCoupling_Free(MRIC);
      return (NULL);
    }
    for (i = 0; i < stages; i++)
    {
      MRIC->group[i] = NULL;
      MRIC->group[i] = (int*)malloc(stages * sizeof(int));
      if (!(MRIC->group[i]))
      {
        MRIStepCoupling_Free(MRIC);
        return (NULL);
      }
      for (j = 0; j < stages; j++) { MRIC->group[i][j] = -1; }
    }
  }

  return (MRIC);
}

/*---------------------------------------------------------------
  Routine to allocate and fill an explicit, implicit, or ImEx
  MRIGARK MRIStepCoupling structure.
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages, int q, int p,
                                       sunrealtype* W, sunrealtype* G,
                                       sunrealtype* c)
{
  int i, j, k;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRIC = NULL;

  /* Check for legal inputs */
  if (nmat < 1 || stages < 1 || !c) { return (NULL); }

  /* Check for method coefficients and set method type */
  if (W && G) { type = MRISTEP_IMEX; }
  else if (W && !G) { type = MRISTEP_EXPLICIT; }
  else if (!W && G) { type = MRISTEP_IMPLICIT; }
  else { return (NULL); }

  /* Allocate MRIStepCoupling structure */
  MRIC = MRIStepCoupling_Alloc(nmat, stages, type);
  if (!MRIC) { return (NULL); }

  /* -------------------------
   * Copy the inputs into MRIC
   * ------------------------- */

  /* Method and embedding order */
  MRIC->q = q;
  MRIC->p = p;

  /* Abscissae */
  for (i = 0; i < stages; i++) { MRIC->c[i] = c[i]; }

  /* Coupling coefficients stored as 1D arrays, based on whether they
     include embedding coefficients */
  if (p == 0)
  {
    /* non-embedded method:  coupling coefficient 1D arrays have
       length nmat * stages * stages, with each stages * stages
       matrix stored in C (row-major) order */
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i < stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX)
          {
            MRIC->W[k][i][j] = W[stages * (stages * k + i) + j];
          }
          if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
          {
            MRIC->G[k][i][j] = G[stages * (stages * k + i) + j];
          }
        }
      }
    }
  }
  else
  {
    /* embedded method:  coupling coefficient 1D arrays have
       length nmat * (stages+1) * stages, with each (stages+1) * stages
       matrix stored in C (row-major) order */
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i <= stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX)
          {
            MRIC->W[k][i][j] = W[(stages + 1) * (stages * k + i) + j];
          }
          if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
          {
            MRIC->G[k][i][j] = G[(stages + 1) * (stages * k + i) + j];
          }
        }
      }
    }
  }
  return (MRIC);
}

/*---------------------------------------------------------------
  Construct the MRIGARK coupling matrix for an MIS method based
  on a given "slow" Butcher table.
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B, int q, int p)
{
  int i, j, stages;
  sunbooleantype padding;
  sunrealtype Asum;
  sunrealtype*** C;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRIC;

  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* Check that input table is non-NULL */
  if (!B) { return (NULL); }

  /* If p>0, check that input table includes embedding coefficients */
  if ((p > 0) && (B->d == NULL)) { return (NULL); }

  /* -----------------------------------
   * Check that the input table is valid
   * ----------------------------------- */

  /* First stage is just old solution */
  Asum = SUNRabs(B->c[0]);
  for (j = 0; j < B->stages; j++) { Asum += SUNRabs(B->A[0][j]); }
  if (Asum > tol) { return (NULL); }

  /* Last stage exceeds 1 */
  if (B->c[B->stages - 1] > ONE + tol) { return (NULL); }

  /* All stages are sorted */
  for (j = 1; j < B->stages; j++)
  {
    if ((B->c[j] - B->c[j - 1]) < -tol) { return (NULL); }
  }

  /* Each stage at most diagonally implicit */
  Asum = ZERO;
  for (i = 0; i < B->stages; i++)
  {
    for (j = i + 1; j < B->stages; j++) { Asum += SUNRabs(B->A[i][j]); }
  }
  if (Asum > tol) { return (NULL); }

  /* -----------------------------------------
   * determine whether the table needs padding
   * ----------------------------------------- */

  padding = SUNFALSE;

  /* Pad if last stage does not equal 1 */
  if (SUNRabs(B->c[B->stages - 1] - ONE) > tol) { padding = SUNTRUE; }

  /* Pad if last row of A does not equal b */
  for (j = 0; j < B->stages; j++)
  {
    if (SUNRabs(B->A[B->stages - 1][j] - B->b[j]) > tol) { padding = SUNTRUE; }
  }

  /* If final stage is implicit and the method contains an embedding,
     we require padding since d != b */
  if ((p > 0) && (SUNRabs(B->A[B->stages - 1][B->stages - 1]) > tol))
  {
    padding = SUNTRUE;
  }
  stages = (padding) ? B->stages + 1 : B->stages;

  /* -------------------------
   * determine the method type
   * ------------------------- */

  /* Check if the table is strictly lower triangular (explicit) */
  type = MRISTEP_EXPLICIT;

  for (i = 0; i < B->stages; i++)
  {
    for (j = i; j < B->stages; j++)
    {
      if (SUNRabs(B->A[i][j]) > tol) { type = MRISTEP_IMPLICIT; }
    }
  }

  /* ----------------------------
   * construct coupling structure
   * ---------------------------- */

  MRIC = MRIStepCoupling_Alloc(1, stages, type);
  if (!MRIC) { return (NULL); }

  /* Copy method/embedding orders */
  MRIC->q = q;
  MRIC->p = p;

  /* Copy abscissae, padding if needed */
  for (i = 0; i < B->stages; i++) { MRIC->c[i] = B->c[i]; }

  if (padding) { MRIC->c[stages - 1] = ONE; }

  /* Construct the coupling table */
  if (type == MRISTEP_EXPLICIT) { C = MRIC->W; }
  else { C = MRIC->G; }

  /* First row is identically zero */
  for (i = 0; i < stages; i++)
  {
    for (j = 0; j < stages; j++) { C[0][i][j] = ZERO; }
  }

  /* Remaining rows = A(2:end,:) - A(1:end-1,:) */
  for (i = 1; i < B->stages; i++)
  {
    for (j = 0; j < B->stages; j++)
    {
      C[0][i][j] = B->A[i][j] - B->A[i - 1][j];
    }
  }

  /* Padded row = b(:) - A(end,:) */
  if (padding)
  {
    for (j = 0; j < B->stages; j++)
    {
      C[0][stages - 1][j] = B->b[j] - B->A[B->stages - 1][j];
    }
  }

  /* Embedded row = d(:) - A(end,:) */
  if (p > 0)
  {
    for (j = 0; j < B->stages; j++)
    {
      C[0][stages][j] = B->d[j] - B->A[B->stages - 1][j];
    }
  }

  return (MRIC);
}

/*---------------------------------------------------------------
  Routine to copy a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC)
{
  int i, j, k, nmat, stages;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRICcopy;

  /* Check for legal input */
  if (!MRIC) { return (NULL); }

  /* Copy method type */
  type = MRIC->type;

  /* Check for stage times */
  if (!(MRIC->c)) { return (NULL); }

  /* Get the number of coupling matrices and stages */
  nmat   = MRIC->nmat;
  stages = MRIC->stages;

  /* Allocate coupling structure */
  MRICcopy = MRIStepCoupling_Alloc(nmat, stages, type);
  if (!MRICcopy) { return (NULL); }

  /* Copy method and embedding orders */
  MRICcopy->q = MRIC->q;
  MRICcopy->p = MRIC->p;

  /* Copy abscissae */
  for (i = 0; i < stages; i++) { MRICcopy->c[i] = MRIC->c[i]; }

  /* Copy explicit coupling matrices W */
  if (MRIC->W)
  {
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i <= stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          MRICcopy->W[k][i][j] = MRIC->W[k][i][j];
        }
      }
    }
  }

  /* Copy implicit coupling matrices G */
  if (MRIC->G)
  {
    for (k = 0; k < nmat; k++)
    {
      for (i = 0; i <= stages; i++)
      {
        for (j = 0; j < stages; j++)
        {
          MRICcopy->G[k][i][j] = MRIC->G[k][i][j];
        }
      }
    }
  }

  /* Copy MERK stage groups */
  if (MRIC->group)
  {
    MRICcopy->ngroup = MRIC->ngroup;
    for (i = 0; i < stages; i++)
    {
      for (j = 0; j < stages; j++)
      {
        MRICcopy->group[i][j] = MRIC->group[i][j];
      }
    }
  }

  return (MRICcopy);
}

/*---------------------------------------------------------------
  Routine to query the MRIStepCoupling structure workspace size
  ---------------------------------------------------------------*/
void MRIStepCoupling_Space(MRIStepCoupling MRIC, sunindextype* liw,
                           sunindextype* lrw)
{
  /* initialize outputs and return if MRIC is not allocated */
  *liw = 0;
  *lrw = 0;
  if (!MRIC) { return; }

  /* fill outputs based on MRIC */
  *liw = 5;
  if (MRIC->c) { *lrw += MRIC->stages; }
  if (MRIC->W) { *lrw += MRIC->nmat * (MRIC->stages + 1) * MRIC->stages; }
  if (MRIC->G) { *lrw += MRIC->nmat * (MRIC->stages + 1) * MRIC->stages; }
  if (MRIC->group) { *liw += MRIC->stages * MRIC->stages; }
}

/*---------------------------------------------------------------
  Routine to free a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Free(MRIStepCoupling MRIC)
{
  int k, i;

  /* Free each field within MRIStepCoupling structure, and then
     free structure itself */
  if (MRIC)
  {
    if (MRIC->c) { free(MRIC->c); }

    if (MRIC->W)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        if (MRIC->W[k])
        {
          for (i = 0; i <= MRIC->stages; i++)
          {
            if (MRIC->W[k][i])
            {
              free(MRIC->W[k][i]);
              MRIC->W[k][i] = NULL;
            }
          }
          free(MRIC->W[k]);
          MRIC->W[k] = NULL;
        }
      }
      free(MRIC->W);
    }

    if (MRIC->G)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        if (MRIC->G[k])
        {
          for (i = 0; i <= MRIC->stages; i++)
          {
            if (MRIC->G[k][i])
            {
              free(MRIC->G[k][i]);
              MRIC->G[k][i] = NULL;
            }
          }
          free(MRIC->G[k]);
          MRIC->G[k] = NULL;
        }
      }
      free(MRIC->G);
    }

    if (MRIC->group)
    {
      for (i = 0; i < MRIC->stages; i++)
      {
        if (MRIC->group[i])
        {
          free(MRIC->group[i]);
          MRIC->group[i] = NULL;
        }
      }
      free(MRIC->group);
    }

    free(MRIC);
  }
}

/*---------------------------------------------------------------
  Routine to print a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Write(MRIStepCoupling MRIC, FILE* outfile)
{
  int i, j, k;

  /* check for valid coupling structure */
  if (!MRIC) { return; }
  if (!(MRIC->W) && !(MRIC->G)) { return; }
  if (!(MRIC->c)) { return; }

  if (MRIC->W)
  {
    for (i = 0; i < MRIC->nmat; i++)
    {
      if (!(MRIC->W[i])) { return; }
      for (j = 0; j <= MRIC->stages; j++)
      {
        if (!(MRIC->W[i][j])) { return; }
      }
    }
  }

  if (MRIC->G)
  {
    for (i = 0; i < MRIC->nmat; i++)
    {
      if (!(MRIC->G[i])) { return; }
      for (j = 0; j <= MRIC->stages; j++)
      {
        if (!(MRIC->G[i][j])) { return; }
      }
    }
  }

  if (MRIC->group)
  {
    for (i = 0; i < MRIC->stages; i++)
    {
      if (!(MRIC->group[i])) { return; }
    }
  }

  switch (MRIC->type)
  {
  case MRISTEP_EXPLICIT: fprintf(outfile, "  type = explicit MRI\n"); break;
  case MRISTEP_IMPLICIT: fprintf(outfile, "  type = implicit MRI\n"); break;
  case MRISTEP_IMEX: fprintf(outfile, "  type = ImEx MRI\n"); break;
  case MRISTEP_MERK: fprintf(outfile, "  type = MERK\n"); break;
  case MRISTEP_SR: fprintf(outfile, "  type = MRISR\n"); break;
  default: fprintf(outfile, "  type = unknown\n");
  }
  fprintf(outfile, "  nmat = %i\n", MRIC->nmat);
  fprintf(outfile, "  stages = %i\n", MRIC->stages);
  fprintf(outfile, "  method order (q) = %i\n", MRIC->q);
  fprintf(outfile, "  embedding order (p) = %i\n", MRIC->p);
  fprintf(outfile, "  c = ");
  for (i = 0; i < MRIC->stages; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "  ", MRIC->c[i]);
  }
  fprintf(outfile, "\n");

  if (MRIC->W)
  {
    for (k = 0; k < MRIC->nmat; k++)
    {
      fprintf(outfile, "  W[%i] = \n", k);
      for (i = 0; i <= MRIC->stages; i++)
      {
        fprintf(outfile, "      ");
        for (j = 0; j < MRIC->stages; j++)
        {
          fprintf(outfile, SUN_FORMAT_E "  ", MRIC->W[k][i][j]);
        }
        fprintf(outfile, "\n");
      }
      fprintf(outfile, "\n");
    }
  }

  if (MRIC->G)
  {
    for (k = 0; k < MRIC->nmat; k++)
    {
      fprintf(outfile, "  G[%i] = \n", k);
      for (i = 0; i <= MRIC->stages; i++)
      {
        fprintf(outfile, "      ");
        for (j = 0; j < MRIC->stages; j++)
        {
          fprintf(outfile, SUN_FORMAT_E "  ", MRIC->G[k][i][j]);
        }
        fprintf(outfile, "\n");
      }
      fprintf(outfile, "\n");
    }
  }

  if (MRIC->group)
  {
    fprintf(outfile, "  ngroup = %i\n", MRIC->ngroup);
    for (i = 0; i < MRIC->ngroup; i++)
    {
      fprintf(outfile, "  group[%i] = ", i);
      for (j = 0; j < MRIC->stages; j++)
      {
        if (MRIC->group[i][j] >= 0)
        {
          fprintf(outfile, "%i ", MRIC->group[i][j]);
        }
      }
      fprintf(outfile, "\n");
    }
  }
}

/* ===========================================================================
 * Private Functions
 * ===========================================================================*/

/* ---------------------------------------------------------------------------
 * Stage type identifier: returns one of the constants
 *
 * MRISTAGE_ERK_FAST    -- standard MIS-like stage
 * MRISTAGE_ERK_NOFAST  -- standard ERK stage
 * MRISTAGE_DIRK_NOFAST -- standard DIRK stage
 * MRISTAGE_DIRK_FAST   -- coupled DIRK + MIS-like stage
 * MRISTAGE_STIFF_ACC   -- "extra" stiffly-accurate stage
 *
 * for each nontrivial stage, or embedding stage, in an MRI-like method.
 * Otherwise (i.e., stage is not in [1,MRIC->stages]), returns
 * ARK_INVALID_TABLE (<0).
 *
 * The stage type is determined by 2 factors (for normal stages):
 * (a) Sum |MRIC->G[:][is][is]| (nonzero => DIRK)
 * (b) MRIC->c[is] - MRIC->c[is-1]  (nonzero => fast)
 * Similar tests are used for embedding stages.
 *
 * Note that MERK and MRI-SR methods do not use the stage-type identifiers,
 * so if those tables are input we just return MRISTAGE_ERK_FAST.
 * ---------------------------------------------------------------------------*/

int mriStepCoupling_GetStageType(MRIStepCoupling MRIC, int is)
{
  int i, j;
  sunbooleantype Gdiag, Grow, Wrow, cdiff;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  if ((is < 0) || (is > MRIC->stages)) { return ARK_INVALID_TABLE; }

  if (is == 0) { return MRISTAGE_FIRST; }

  /* report MRISTAGE_ERK_FAST for MERK and MRI-SR methods */
  if ((MRIC->type == MRISTEP_SR) || (MRIC->type == MRISTEP_MERK))
  {
    return (MRISTAGE_ERK_FAST);
  }

  /* separately handle an embedding "stage" from normal stages */
  if (is < MRIC->stages)
  { /* normal */
    Gdiag = Grow = Wrow = cdiff = SUNFALSE;
    if (MRIC->G)
    {
      for (i = 0; i < MRIC->nmat; i++)
      {
        Gdiag = Gdiag || (SUNRabs(MRIC->G[i][is][is]) > tol);
        for (j = 0; j < MRIC->stages; j++)
        {
          Grow = Grow || (SUNRabs(MRIC->G[i][is][j]) > tol);
        }
      }
    }
    if (MRIC->W)
    {
      for (i = 0; i < MRIC->nmat; i++)
      {
        for (j = 0; j < MRIC->stages; j++)
        {
          Wrow = Wrow || (SUNRabs(MRIC->W[i][is][j]) > tol);
        }
      }
    }

    /* abscissae difference */
    cdiff = (SUNRabs(MRIC->c[is] - MRIC->c[is - 1]) > tol);
  }
  else
  { /* embedding */
    Gdiag = Grow = Wrow = cdiff = SUNFALSE;
    if (MRIC->G)
    {
      for (i = 0; i < MRIC->nmat; i++)
      {
        Gdiag = Gdiag || (SUNRabs(MRIC->G[i][is][is - 1]) > tol);
        for (j = 0; j < MRIC->stages; j++)
        {
          Grow = Grow || (SUNRabs(MRIC->G[i][is][j]) > tol);
        }
      }
    }
    if (MRIC->W)
    {
      for (i = 0; i < MRIC->nmat; i++)
      {
        for (j = 0; j < MRIC->stages; j++)
        {
          Wrow = Wrow || (SUNRabs(MRIC->W[i][is][j]) > tol);
        }
      }
    }
    cdiff = (SUNRabs(MRIC->c[is - 1] - MRIC->c[is - 2]) > tol);
  }

  /* make determination */
  if (!(Gdiag || Grow || Wrow || cdiff) && (is > 0))
  { /* stiffly-accurate stage */
    return (MRISTAGE_STIFF_ACC);
  }
  if (Gdiag)
  { /* DIRK */
    if (cdiff)
    { /* Fast */
      return (MRISTAGE_DIRK_FAST);
    }
    else { return (MRISTAGE_DIRK_NOFAST); }
  }
  else
  { /* ERK */
    if (cdiff)
    { /* Fast */
      return (MRISTAGE_ERK_FAST);
    }
    else { return (MRISTAGE_ERK_NOFAST); }
  }
}

/* ---------------------------------------------------------------------------
 * Computes the stage RHS vector storage maps. With repeated abscissae the
 * first stage of the pair generally corresponds to a column of zeros and so
 * does not need to be computed and stored. The stage_map indicates if the RHS
 * needs to be computed and where to store it i.e., stage_map[i] > -1.
 *
 * Note: for MERK and MRI-SR methods, this should be an "identity" map, and all
 * stage vectors should be allocated.
 * ---------------------------------------------------------------------------*/

int mriStepCoupling_GetStageMap(MRIStepCoupling MRIC, int* stage_map,
                                int* nstages_active)
{
  int i, j, k, idx;
  sunrealtype Wsum, Gsum;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* ----------------------
   * Check for valid inputs
   * ---------------------- */

  if (!MRIC) { return (ARK_ILL_INPUT); }
  if (!(MRIC->W) && !(MRIC->G)) { return (ARK_ILL_INPUT); }
  if (!stage_map || !nstages_active) { return (ARK_ILL_INPUT); }

  /* -------------------------------------------
   * MERK and MRI-SR have "identity" storage map
   * ------------------------------------------- */

  if ((MRIC->type == MRISTEP_MERK) || (MRIC->type == MRISTEP_SR))
  {
    /* Number of stage RHS vectors active */
    *nstages_active = MRIC->stages;

    /* Create an identity map (all columns are non-zero) */
    for (j = 0; j < MRIC->stages; j++) { stage_map[j] = j; }
    return (ARK_SUCCESS);
  }

  /* ----------------------------------------
   * Compute storage map for MRI-GARK methods
   * ---------------------------------------- */

  /* Number of stage RHS vectors active */
  *nstages_active = 0;

  /* Initial storage index */
  idx = 0;

  /* Check if a stage corresponds to a column of zeros for all coupling
   * matrices by computing the column sums */
  for (j = 0; j < MRIC->stages; j++)
  {
    Wsum = ZERO;
    Gsum = ZERO;

    if (MRIC->W)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        for (i = 0; i <= MRIC->stages; i++)
        {
          Wsum += SUNRabs(MRIC->W[k][i][j]);
        }
      }
    }

    if (MRIC->G)
    {
      for (k = 0; k < MRIC->nmat; k++)
      {
        for (i = 0; i <= MRIC->stages; i++)
        {
          Gsum += SUNRabs(MRIC->G[k][i][j]);
        }
      }
    }

    if (Wsum > tol || Gsum > tol)
    {
      stage_map[j] = idx;
      idx++;
    }
    else { stage_map[j] = -1; }
  }

  /* Check and set number of stage RHS vectors active */
  if (idx < 1) { return (ARK_ILL_INPUT); }

  *nstages_active = idx;

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
