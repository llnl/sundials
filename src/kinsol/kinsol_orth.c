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
 * Orthogonalization utilities
 * ---------------------------------------------------------------------------*/

#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

#include "kinsol_impl.h"

int KINInitOrth(KINMem kin_mem)
{
  // Do we need to (re)allocate the orthogonalization workspace?
  sunbooleantype allocate = kin_mem->kin_m_aa > kin_mem->kin_orth_aa_alloc;

  if (allocate)
  {
    // Free any existing workspace allocations
    KINFreeOrth(kin_mem);

    // Template vector for creating clones
    N_Vector tmpl = kin_mem->kin_unew;

    // Update the AA workspace size
    kin_mem->kin_orth_aa_alloc = kin_mem->kin_m_aa;

    // Structure of orthogonalization data for QR solve
    kin_mem->kin_qr_data = (SUNQRData)malloc(sizeof(*kin_mem->kin_qr_data));
    if (kin_mem->kin_qr_data == NULL)
    {
      KINFreeOrth(kin_mem);
      KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
      return KIN_MEM_FAIL;
    }

    if (kin_mem->kin_orth_aa != KIN_ORTH_MGS)
    {
      kin_mem->kin_vtemp3 = N_VClone(tmpl); // Orth owns
      if (kin_mem->kin_vtemp3 == NULL)
      {
        KINFreeOrth(kin_mem);
        KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
        return KIN_MEM_FAIL;
      }
    }

    if (kin_mem->kin_orth_aa == KIN_ORTH_ICWY)
    {
      // T matrix for ICWY
      kin_mem->kin_T_aa = (sunrealtype*)malloc(
        ((kin_mem->kin_m_aa * kin_mem->kin_m_aa)) * sizeof(sunrealtype));
      if (kin_mem->kin_T_aa == NULL)
      {
        KINFreeOrth(kin_mem);
        KINProcessError(kin_mem, 0, __LINE__, __func__, __FILE__, MSG_MEM_FAIL);
        return KIN_MEM_FAIL;
      }
    }
  }

  // Does the vector support dot product with single buffer reductions
  kin_mem->kin_dot_prod_sb = SUNFALSE;
  if ((kin_mem->kin_unew->ops->nvdotprodlocal ||
       kin_mem->kin_unew->ops->nvdotprodmultilocal) &&
      kin_mem->kin_unew->ops->nvdotprodmultiallreduce)
  {
    kin_mem->kin_dot_prod_sb = SUNTRUE;
  }

  // Initialize the QRData and set the QRAdd function
  if (kin_mem->kin_orth_aa == KIN_ORTH_MGS)
  {
    kin_mem->kin_qr_func        = (SUNQRAddFn)SUNQRAdd_MGS;
    kin_mem->kin_qr_data->vtemp = kin_mem->kin_vtemp2; // KINSOL owns
  }
  else if (kin_mem->kin_orth_aa == KIN_ORTH_ICWY)
  {
    if (kin_mem->kin_dot_prod_sb)
    {
      kin_mem->kin_qr_func = (SUNQRAddFn)SUNQRAdd_ICWY_SB;
    }
    else { kin_mem->kin_qr_func = (SUNQRAddFn)SUNQRAdd_ICWY; }
    kin_mem->kin_qr_data->vtemp      = kin_mem->kin_vtemp2; // KINSOL owns
    kin_mem->kin_qr_data->vtemp2     = kin_mem->kin_vtemp3; // Orth owns
    kin_mem->kin_qr_data->temp_array = kin_mem->kin_T_aa;   // Orth owns
  }
  else if (kin_mem->kin_orth_aa == KIN_ORTH_CGS2)
  {
    kin_mem->kin_qr_func             = (SUNQRAddFn)SUNQRAdd_CGS2;
    kin_mem->kin_qr_data->vtemp      = kin_mem->kin_vtemp2; // KINSOL owns
    kin_mem->kin_qr_data->vtemp2     = kin_mem->kin_vtemp3; // Orth owns
    kin_mem->kin_qr_data->temp_array = kin_mem->kin_cv;     // AA owns
  }
  else if (kin_mem->kin_orth_aa == KIN_ORTH_DCGS2)
  {
    if (kin_mem->kin_dot_prod_sb)
    {
      kin_mem->kin_qr_func = (SUNQRAddFn)SUNQRAdd_DCGS2_SB;
    }
    else { kin_mem->kin_qr_func = (SUNQRAddFn)SUNQRAdd_DCGS2; }
    kin_mem->kin_qr_data->vtemp      = kin_mem->kin_vtemp2; // KINSOL owns
    kin_mem->kin_qr_data->vtemp2     = kin_mem->kin_vtemp3; // Orth owns
    kin_mem->kin_qr_data->temp_array = kin_mem->kin_cv;     // AA owns
  }

  return KIN_SUCCESS;
}

void KINFreeOrth(KINMem kin_mem)
{
  if (kin_mem->kin_qr_data)
  {
    free(kin_mem->kin_qr_data);
    kin_mem->kin_qr_data = NULL;
  }

  if (kin_mem->kin_vtemp3)
  {
    N_VDestroy(kin_mem->kin_vtemp3);
    kin_mem->kin_vtemp3 = NULL;
  }

  if (kin_mem->kin_T_aa)
  {
    free(kin_mem->kin_T_aa);
    kin_mem->kin_T_aa = NULL;
  }

  // Reset AA workspace size
  kin_mem->kin_orth_aa_alloc = 0;
}
