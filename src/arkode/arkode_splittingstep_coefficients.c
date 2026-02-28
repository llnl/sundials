/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This is the implementation file for splitting coefficients
 *--------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include <arkode/arkode_splittingstep.h>

#include "arkode_impl.h"

/*---------------------------------------------------------------
  Routine to allocate splitting coefficients with zero values for
  alpha and beta
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_Alloc(
  const int sequential_methods, const int stages, const int partitions)
{
  if (sequential_methods < 1 || stages < 1 || partitions < 1) { return NULL; }

  SplittingStepCoefficients coefficients =
    (SplittingStepCoefficients)malloc(sizeof(*coefficients));
  if (coefficients == NULL) { return NULL; }

  coefficients->sequential_methods = sequential_methods;
  coefficients->stages             = stages;
  coefficients->partitions         = partitions;
  coefficients->order              = 0;

  coefficients->alpha = (sunrealtype*)calloc(sequential_methods,
                                             sizeof(*coefficients->alpha));
  if (coefficients->alpha == NULL)
  {
    SplittingStepCoefficients_Destroy(&coefficients);
    return NULL;
  }

  /* In order for beta to be indexed like beta[i][j][k], we can't use a single
     malloc. The j index requires allocating a matrix of pointers into beta. The
     i index requires allocating an array of pointers into that matrix. */

  /* Array of pointers for index i */
  coefficients->beta =
    (sunrealtype***)malloc(sequential_methods * sizeof(*coefficients->beta));
  if (coefficients->beta == NULL)
  {
    SplittingStepCoefficients_Destroy(&coefficients);
    return NULL;
  }

  /* Matrix of pointers for index j */
  sunrealtype** beta_cols = (sunrealtype**)malloc(
    sequential_methods * (stages + 1) * sizeof(*beta_cols));
  if (beta_cols == NULL)
  {
    SplittingStepCoefficients_Destroy(&coefficients);
    return NULL;
  }

  for (int i = 0; i < sequential_methods; i++)
  {
    coefficients->beta[i] = &beta_cols[i * (stages + 1)];
  }

  /* Contiguous memory to store the beta tensor -- use calloc so only non-zero
   * coefficients need to be set */
  sunrealtype* beta_mem =
    (sunrealtype*)calloc(sequential_methods * (stages + 1) * partitions,
                         sizeof(*beta_mem));
  if (beta_mem == NULL)
  {
    SplittingStepCoefficients_Destroy(&coefficients);
    return NULL;
  }

  /* Set pointers into the beta tensor */
  for (int i = 0; i < sequential_methods; i++)
  {
    for (int j = 0; j <= stages; j++)
    {
      coefficients->beta[i][j] = &beta_mem[(i * (stages + 1) + j) * partitions];
    }
  }

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to create splitting coefficients which performs a copy
  of the alpha and beta parameters
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_Create(
  const int sequential_methods, const int stages, const int partitions,
  const int order, sunrealtype* const alpha, sunrealtype* const beta)
{
  if (alpha == NULL || beta == NULL || order < 1) { return NULL; }

  SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Alloc(sequential_methods, stages, partitions);
  if (coefficients == NULL) { return NULL; }

  coefficients->order = order;
  memcpy(coefficients->alpha, alpha, sequential_methods * sizeof(sunrealtype));
  memcpy(coefficients->beta[0][0], beta,
         sequential_methods * (stages + 1) * partitions * sizeof(sunrealtype));

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to free splitting coefficients
  ---------------------------------------------------------------*/
void SplittingStepCoefficients_Destroy(SplittingStepCoefficients* coefficients)
{
  if (coefficients == NULL || *coefficients == NULL) { return; }

  SplittingStepCoefficients coeffs = *coefficients;
  if (coeffs->alpha != NULL) { free(coeffs->alpha); }
  if (coeffs->beta != NULL)
  {
    if (coeffs->beta[0] != NULL)
    {
      if (coeffs->beta[0][0] != NULL) { free(coeffs->beta[0][0]); }
      free(coeffs->beta[0]);
    }
    free(coeffs->beta);
  }
  free(coeffs);
  *coefficients = NULL;
}

/*---------------------------------------------------------------
  Routine to create a copy of splitting coefficients
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_Copy(
  const SplittingStepCoefficients coefficients)
{
  if (coefficients == NULL) { return NULL; }

  SplittingStepCoefficients coefficientsCopy =
    SplittingStepCoefficients_Alloc(coefficients->sequential_methods,
                                    coefficients->stages,
                                    coefficients->partitions);
  if (coefficientsCopy == NULL) { return NULL; }

  coefficientsCopy->order = coefficients->order;
  memcpy(coefficientsCopy->alpha, coefficients->alpha,
         coefficients->sequential_methods * sizeof(sunrealtype));

  /* beta[0][0] points to the contiguous memory allocation, so we can copy it
     with a single memcpy */
  memcpy(coefficientsCopy->beta[0][0], coefficients->beta[0][0],
         coefficients->sequential_methods * (coefficients->stages + 1) *
           coefficients->partitions * sizeof(sunrealtype));

  return coefficientsCopy;
}

/*---------------------------------------------------------------
  Routine to load coefficients from an ID
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_LoadCoefficients(
  const ARKODE_SplittingCoefficientsID method)
{
  switch (method)
  {
#define ARK_SPLITTING_COEFFICIENTS(name, coeff) \
  case name: coeff break;
#include "arkode_splittingstep_coefficients.def"
#undef ARK_SPLITTING_COEFFICIENTS

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown splitting coefficients");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Routine to load coefficients using a string representation of
  an enum entry in ARKODE_SplittingCoefficientsID
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_LoadCoefficientsByName(
  const char* const method)
{
#define ARK_SPLITTING_COEFFICIENTS(name, coeff) \
  if (strcmp(#name, method) == 0) coeff
#include "arkode_splittingstep_coefficients.def"
#undef ARK_SPLITTING_COEFFICIENTS

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown splitting coefficients");

  return NULL;
}

/*---------------------------------------------------------------
  Routine to convert a coefficient enum value to its string
  representation
  ---------------------------------------------------------------*/
const char* SplittingStepCoefficients_IDToName(
  const ARKODE_SplittingCoefficientsID id)
{
  /* Use X-macro to test each coefficient name */
  switch (id)
  {
#define ARK_SPLITTING_COEFFICIENTS(name, coeff) \
  case name: return #name;
#include "arkode_splittingstep_coefficients.def"
#undef ARK_SPLITTING_COEFFICIENTS

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown splitting coefficients");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Routine to construct the standard Lie-Trotter splitting
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_LieTrotter(const int partitions)
{
  const SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Alloc(1, 1, partitions);
  if (coefficients == NULL) { return NULL; }

  coefficients->order    = 1;
  coefficients->alpha[0] = SUN_RCONST(1.0);
  for (int i = 0; i < partitions; i++)
  {
    coefficients->beta[0][1][i] = SUN_RCONST(1.0);
  }

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to construct the standard Stang splitting
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_Strang(const int partitions)
{
  return SplittingStepCoefficients_TripleJump(partitions, 2);
}

/*---------------------------------------------------------------
  Routine to construct a parallel splitting method
  Phi_1(h) + Phi_2(h) + ... + Phi_p(h) - (p - 1) * y_n
  where Phi_i is the flow of partition i and p = partitions.
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_Parallel(const int partitions)
{
  const SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Alloc(partitions + 1, 1, partitions);
  if (coefficients == NULL) { return NULL; }

  coefficients->order = 1;
  for (int i = 0; i < partitions; i++)
  {
    coefficients->alpha[i]      = SUN_RCONST(1.0);
    coefficients->beta[i][1][i] = SUN_RCONST(1.0);
  }

  coefficients->alpha[partitions] = 1 - partitions;

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to construct a symmetric parallel splitting which is
  the average of the Lie-Trotter method and its adjoint
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_SymmetricParallel(
  const int partitions)
{
  const SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Alloc(2, partitions, partitions);
  if (coefficients == NULL) { return NULL; }

  coefficients->order    = 2;
  coefficients->alpha[0] = SUN_RCONST(0.5);
  coefficients->alpha[1] = SUN_RCONST(0.5);

  for (int i = 0; i < partitions; i++)
  {
    coefficients->beta[0][partitions][i] = SUN_RCONST(1.0);
    for (int j = partitions - i - 1; j < partitions; j++)
    {
      coefficients->beta[1][i + 1][j] = SUN_RCONST(1.0);
    }
  }

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to construct a 3rd order method of Suzuki of the form
  L(p1 h) * L*(p2 h) * L(p3 h) * L*(p4 h) * L(p5 h)
  where L is a Lie-Trotter splitting and L* is its adjoint.
  Composition is denoted by *.
  ---------------------------------------------------------------*/
SplittingStepCoefficients SplittingStepCoefficients_ThirdOrderSuzuki(
  const int partitions)
{
  const SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Alloc(1, 2 * partitions - 1, partitions);
  if (coefficients == NULL) { return NULL; }

  coefficients->order    = 3;
  coefficients->alpha[0] = SUN_RCONST(1.0);

  for (int i = 1; i < partitions; i++)
  {
    for (int j = 0; j < partitions; j++)
    {
      // Constants from https://doi.org/10.1143/JPSJ.61.3015 pg. 3019
      const sunrealtype p1 =
        SUN_RCONST(0.2683300957817599249569552299254991394812);
      const sunrealtype p2 =
        SUN_RCONST(0.6513314272356399320939424082278836500821);

      coefficients->beta[0][i][j] = i + j < partitions ? p1 : (p1 + p2);
      coefficients->beta[0][partitions + i - 1][j] =
        SUN_RCONST(1.0) - (i + j < partitions ? (p1 + p2) : p1);
    }
  }

  for (int i = 0; i < partitions; i++)
  {
    coefficients->beta[0][2 * partitions - 1][i] = SUN_RCONST(1.0);
  }

  return coefficients;
}

/*---------------------------------------------------------------
  Routine to construct a composition method of the form
  S(gamma_0 h)^c * S(gamma_1 h) * S(gamma_0)^c
  where S is a lower order splitting (with Stang as the base case),
  * and ^ denote composition, and c = composition_stages. This
  covers both the triple jump (c=1) and Suzuki fractal (c=2).
  ---------------------------------------------------------------*/
static sunrealtype* const* SplittingStepCoefficients_ComposeStrangHelper(
  const int partitions, const int order, const int composition_stages,
  const sunrealtype start, const sunrealtype end, sunrealtype* const* const beta)
{
  const sunrealtype diff = end - start;
  if (order == 2)
  {
    /* The base case is an order 2 Strang splitting */
    const sunrealtype mid = start + diff / SUN_RCONST(2.0);
    for (int j = 1; j <= partitions; j++)
    {
      for (int k = 0; k < partitions; k++)
      {
        beta[j][k] = (k + j < partitions) ? mid : end;
      }
    }

    return &beta[partitions - 1];
  }

  sunrealtype* const* beta_cur = beta;
  sunrealtype start_cur        = start;
  /* This is essentially the gamma coefficient from Geometric Numerical
   * Integration (https://doi.org/10.1007/3-540-30666-8) pg 44-45 scaled by the
   * current interval */
  const sunrealtype gamma =
    diff / (composition_stages - 1 -
            SUNRpowerR(composition_stages - 1, SUN_RCONST(1.0) / (order - 1)));
  for (int i = 1; i <= composition_stages; i++)
  {
    /* To avoid roundoff issues, this ensures end_cur=1 for the last value of i*/
    const sunrealtype end_cur = 2 * i < composition_stages
                                  ? (start + i * gamma)
                                  : (end + (i - composition_stages) * gamma);
    /* Recursively generate coefficients and shift beta_cur */
    beta_cur  = SplittingStepCoefficients_ComposeStrangHelper(partitions,
                                                              order - 2,
                                                              composition_stages,
                                                              start_cur, end_cur,
                                                              beta_cur);
    start_cur = end_cur;
  }

  return beta_cur;
}

/*---------------------------------------------------------------
  Routine which does validation and setup before calling
  SplittingStepCoefficients_ComposeStrangHelper to fill in the
  beta coefficients
  ---------------------------------------------------------------*/
static SplittingStepCoefficients SplittingStepCoefficients_ComposeStrang(
  const int partitions, const int order, const int composition_stages)
{
  if (order < 2 || order % 2 != 0)
  {
    // Only even orders allowed
    return NULL;
  }

  const int stages = 1 + (partitions - 1) *
                           SUNIpowerI(composition_stages, order / 2 - 1);
  const SplittingStepCoefficients coefficients =
    SplittingStepCoefficients_Alloc(1, stages, partitions);
  if (coefficients == NULL) { return NULL; }

  coefficients->order    = order;
  coefficients->alpha[0] = SUN_RCONST(1.0);

  SplittingStepCoefficients_ComposeStrangHelper(partitions, order,
                                                composition_stages,
                                                SUN_RCONST(0.0), SUN_RCONST(1.0),
                                                coefficients->beta[0]);

  return coefficients;
}

SplittingStepCoefficients SplittingStepCoefficients_TripleJump(const int partitions,
                                                               const int order)
{
  return SplittingStepCoefficients_ComposeStrang(partitions, order, 3);
}

SplittingStepCoefficients SplittingStepCoefficients_SuzukiFractal(
  const int partitions, const int order)
{
  return SplittingStepCoefficients_ComposeStrang(partitions, order, 5);
}

/*---------------------------------------------------------------
  Routine to print a splitting coefficient structure
  ---------------------------------------------------------------*/
void SplittingStepCoefficients_Write(const SplittingStepCoefficients coefficients,
                                     FILE* const outfile)
{
  if (outfile == NULL || coefficients == NULL || coefficients->alpha == NULL ||
      coefficients->beta == NULL || coefficients->beta[0] == NULL ||
      coefficients->beta[0][0] == NULL)
  {
    return;
  }

  fprintf(outfile, "  sequential methods = %i\n",
          coefficients->sequential_methods);
  fprintf(outfile, "  stages = %i\n", coefficients->stages);
  fprintf(outfile, "  partitions = %i\n", coefficients->partitions);
  fprintf(outfile, "  order = %i\n", coefficients->order);
  fprintf(outfile, "  alpha = ");
  for (int i = 0; i < coefficients->sequential_methods; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "  ", coefficients->alpha[i]);
  }
  fprintf(outfile, "\n");

  for (int i = 0; i < coefficients->sequential_methods; i++)
  {
    fprintf(outfile, "  beta[%i] = \n", i);
    for (int j = 0; j <= coefficients->stages; j++)
    {
      fprintf(outfile, "      ");
      for (int k = 0; k < coefficients->partitions; k++)
      {
        fprintf(outfile, SUN_FORMAT_E "  ", coefficients->beta[i][j][k]);
      }
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }
}
