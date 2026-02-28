/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Anderson acceleration utilities
 * ---------------------------------------------------------------------------*/

#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

#include "kinsol_impl.h"

int KINInitAA(KINMem kin_mem)
{
  // Limit the acceleration space size
  if (kin_mem->kin_m_aa >= kin_mem->kin_mxiter)
  {
    kin_mem->kin_m_aa = kin_mem->kin_mxiter - 1;
  }

  // Initialize the current depth
  kin_mem->kin_current_depth = 0;

  // Do we need to (re)allocate the AA workspace?
  sunbooleantype allocate = kin_mem->kin_m_aa > kin_mem->kin_m_aa_alloc;

  if (allocate)
  {
    // Free any existing workspace allocations
    KINFreeAA(kin_mem);

    // Template vector for creating clones
    N_Vector tmpl = kin_mem->kin_unew;

    // Update the AA workspace size
    kin_mem->kin_m_aa_alloc = kin_mem->kin_m_aa;

    // Array of acceleration weights
    kin_mem->kin_gamma_aa =
      (sunrealtype*)malloc(kin_mem->kin_m_aa * sizeof(sunrealtype));
    if (kin_mem->kin_gamma_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // R matrix for QR factorization
    kin_mem->kin_R_aa = (sunrealtype*)malloc(
      (kin_mem->kin_m_aa * kin_mem->kin_m_aa) * sizeof(sunrealtype));
    if (kin_mem->kin_R_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // Q matrix for QR factorization
    kin_mem->kin_q_aa = N_VCloneVectorArray((int)kin_mem->kin_m_aa, tmpl);
    if (kin_mem->kin_q_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // History of residual vector differences
    kin_mem->kin_df_aa = N_VCloneVectorArray((int)kin_mem->kin_m_aa, tmpl);
    if (kin_mem->kin_df_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // History of fixed point function vector differences
    kin_mem->kin_dg_aa = N_VCloneVectorArray((int)kin_mem->kin_m_aa, tmpl);
    if (kin_mem->kin_dg_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // Previous residual vector, F(u_{i-1}) = G(u_{i-1}) - u_{i-1}
    kin_mem->kin_fold_aa = N_VClone(tmpl);
    if (kin_mem->kin_fold_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // Previous fixed point function vector, G(u_{i-1})
    kin_mem->kin_gold_aa = N_VClone(tmpl);
    if (kin_mem->kin_gold_aa == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // Workspace array for fused operation constants
    kin_mem->kin_cv =
      (sunrealtype*)malloc(2 * (kin_mem->kin_m_aa + 1) * sizeof(sunrealtype));
    if (kin_mem->kin_cv == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    // Workspace array for fused operation vectors
    kin_mem->kin_Xv =
      (N_Vector*)malloc(2 * (kin_mem->kin_m_aa + 1) * sizeof(N_Vector));
    if (kin_mem->kin_Xv == NULL)
    {
      KINFreeAA(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }
  }

  return KIN_SUCCESS;
}

void KINFreeAA(KINMem kin_mem)
{
  if (kin_mem->kin_gamma_aa)
  {
    free(kin_mem->kin_gamma_aa);
    kin_mem->kin_gamma_aa = NULL;
  }

  if (kin_mem->kin_R_aa)
  {
    free(kin_mem->kin_R_aa);
    kin_mem->kin_R_aa = NULL;
  }

  if (kin_mem->kin_q_aa)
  {
    N_VDestroyVectorArray(kin_mem->kin_q_aa, (int)kin_mem->kin_m_aa);
    kin_mem->kin_q_aa = NULL;
  }

  if (kin_mem->kin_df_aa)
  {
    N_VDestroyVectorArray(kin_mem->kin_df_aa, (int)kin_mem->kin_m_aa);
    kin_mem->kin_df_aa = NULL;
  }

  if (kin_mem->kin_dg_aa)
  {
    N_VDestroyVectorArray(kin_mem->kin_dg_aa, (int)kin_mem->kin_m_aa);
    kin_mem->kin_dg_aa = NULL;
  }

  if (kin_mem->kin_fold_aa)
  {
    N_VDestroy(kin_mem->kin_fold_aa);
    kin_mem->kin_fold_aa = NULL;
  }

  if (kin_mem->kin_gold_aa)
  {
    N_VDestroy(kin_mem->kin_gold_aa);
    kin_mem->kin_gold_aa = NULL;
  }

  if (kin_mem->kin_cv)
  {
    free(kin_mem->kin_cv);
    kin_mem->kin_cv = NULL;
  }

  if (kin_mem->kin_Xv)
  {
    free(kin_mem->kin_Xv);
    kin_mem->kin_Xv = NULL;
  }

  // Reset AA workspace size
  kin_mem->kin_m_aa_alloc = 0;
}
