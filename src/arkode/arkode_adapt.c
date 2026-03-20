/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the implementation file for ARKODE's time step
 * adaptivity utilities.
 *--------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"

/*---------------------------------------------------------------
  arkAdaptInit:

  This routine creates and sets default values in an
  ARKodeHAdaptMem structure.  This returns a non-NULL structure
  if no errors occurred, or a NULL value otherwise.
  ---------------------------------------------------------------*/
ARKodeHAdaptMem arkAdaptInit(void)
{
  ARKodeHAdaptMem hadapt_mem;

  /* allocate structure */
  hadapt_mem = (ARKodeHAdaptMem)malloc(sizeof(struct ARKodeHAdaptMemRec));
  if (hadapt_mem == NULL) { return (NULL); }

  /* initialize values (default parameters are set in ARKodeSetDefaults) */
  memset(hadapt_mem, 0, sizeof(struct ARKodeHAdaptMemRec));
  hadapt_mem->nst_acc = 0;
  hadapt_mem->nst_exp = 0;
  return (hadapt_mem);
}

/*---------------------------------------------------------------
  arkPrintAdaptMem

  This routine outputs the time step adaptivity memory structure
  to a specified file pointer.
  ---------------------------------------------------------------*/
void arkPrintAdaptMem(ARKodeHAdaptMem hadapt_mem, FILE* outfile)
{
  if (hadapt_mem != NULL)
  {
    fprintf(outfile, "ark_hadapt: etamax = " SUN_FORMAT_G "\n",
            hadapt_mem->etamax);
    fprintf(outfile, "ark_hadapt: etamx1 = " SUN_FORMAT_G "\n",
            hadapt_mem->etamx1);
    fprintf(outfile, "ark_hadapt: etamxf = " SUN_FORMAT_G "\n",
            hadapt_mem->etamxf);
    fprintf(outfile, "ark_hadapt: etamin = " SUN_FORMAT_G "\n",
            hadapt_mem->etamin);
    fprintf(outfile, "ark_hadapt: small_nef = %i\n", hadapt_mem->small_nef);
    fprintf(outfile, "ark_hadapt: etacf = " SUN_FORMAT_G "\n", hadapt_mem->etacf);
    fprintf(outfile, "ark_hadapt: cfl = " SUN_FORMAT_G "\n", hadapt_mem->cfl);
    fprintf(outfile, "ark_hadapt: safety = " SUN_FORMAT_G "\n",
            hadapt_mem->safety);
    fprintf(outfile, "ark_hadapt: growth = " SUN_FORMAT_G "\n",
            hadapt_mem->growth);
    fprintf(outfile, "ark_hadapt: lbound = " SUN_FORMAT_G "\n",
            hadapt_mem->lbound);
    fprintf(outfile, "ark_hadapt: ubound = " SUN_FORMAT_G "\n",
            hadapt_mem->ubound);
    fprintf(outfile, "ark_hadapt: nst_acc = %li\n", hadapt_mem->nst_acc);
    fprintf(outfile, "ark_hadapt: nst_exp = %li\n", hadapt_mem->nst_exp);
    fprintf(outfile, "ark_hadapt: pq = %i\n", hadapt_mem->pq);
    fprintf(outfile, "ark_hadapt: p = %i\n", hadapt_mem->p);
    fprintf(outfile, "ark_hadapt: q = %i\n", hadapt_mem->q);
    fprintf(outfile, "ark_hadapt: adjust = %i\n", hadapt_mem->adjust);
    if (hadapt_mem->expstab == NULL)
    {
      fprintf(outfile,
              "  ark_hadapt: No explicit stability function supplied\n");
    }
    else
    {
      fprintf(outfile,
              "  ark_hadapt: User provided explicit stability function\n");
      fprintf(outfile, "  ark_hadapt: stability function data pointer = %p\n",
              hadapt_mem->estab_data);
    }
    if (hadapt_mem->hcontroller != NULL)
    {
      (void)SUNAdaptController_Write(hadapt_mem->hcontroller, outfile);
    }
  }
}

/*---------------------------------------------------------------
  arkAdapt is the time step adaptivity wrapper function.  This
  computes and sets the value of ark_eta inside of the ARKodeMem
  data structure.
  ---------------------------------------------------------------*/
int arkAdapt(ARKodeMem ark_mem, ARKodeHAdaptMem hadapt_mem, N_Vector ycur,
             sunrealtype tcur, sunrealtype hcur, sunrealtype dsm)
{
  int retval;
  sunrealtype h_acc;
  int controller_order;

  /* Return with no stepsize adjustment if the controller is NULL */
  if (hadapt_mem->hcontroller == NULL)
  {
    ark_mem->eta = ONE;
    return (ARK_SUCCESS);
  }

  /* Request error-based step size from adaptivity controller */
  if (hadapt_mem->pq == 0)
  {
    controller_order = hadapt_mem->p + hadapt_mem->adjust;
  }
  else if (hadapt_mem->pq == 1)
  {
    controller_order = hadapt_mem->q + hadapt_mem->adjust;
  }
  else
  {
    controller_order = SUNMIN(hadapt_mem->p, hadapt_mem->q) + hadapt_mem->adjust;
  }
  retval = SUNAdaptController_EstimateStep(hadapt_mem->hcontroller, hcur,
                                           controller_order, dsm, &h_acc);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__, __FILE__,
                    "SUNAdaptController_EstimateStep failure.");
    return (ARK_CONTROLLER_ERR);
  }

  SUNLogDebug(ARK_LOGGER, "new-step-before-bounds", "h_acc = " SUN_FORMAT_G,
              h_acc);

  /* enforce safety factors */
  h_acc *= hadapt_mem->safety;

  /* enforce maximum bound on time step growth */
  h_acc = SUNMIN(SUNRabs(h_acc), SUNRabs(hadapt_mem->etamax * hcur));

  /* enforce minimum bound time step reduction */
  h_acc = SUNMAX(h_acc, SUNRabs(hadapt_mem->etamin * hcur));

  if (hadapt_mem->expstab != NULL)
  {
    sunrealtype h_cfl = ZERO;
    retval = hadapt_mem->expstab(ycur, tcur, &h_cfl, hadapt_mem->estab_data);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in explicit stability function.");
      return (ARK_ILL_INPUT);
    }

    h_cfl *= hadapt_mem->cfl;
    SUNLogDebug(ARK_LOGGER, "new-step-cfl", "h_cfl = " SUN_FORMAT_G, h_cfl);

    if (h_cfl > ZERO && h_cfl < h_acc)
    {
      hadapt_mem->nst_exp++;
      h_acc = h_cfl;
    }
    else { hadapt_mem->nst_acc++; }
  }
  else { hadapt_mem->nst_acc++; }

  /* enforce adaptivity bounds to retain Jacobian/preconditioner accuracy */
  if (dsm <= ONE)
  {
    if ((h_acc > SUNRabs(hcur * hadapt_mem->lbound * ONEMSM)) &&
        (h_acc < SUNRabs(hcur * hadapt_mem->ubound * ONEPSM)))
    {
      h_acc = hcur;
    }
  }
  h_acc = SUNRcopysign(h_acc, hcur);

  SUNLogDebug(ARK_LOGGER, "new-step-after-max-min-bounds",
              "h_acc = " SUN_FORMAT_G, h_acc);

  /* set basic value of ark_eta */
  ark_mem->eta = h_acc / hcur;

  /* enforce minimum time step size */
  ark_mem->eta = SUNMAX(ark_mem->eta, ark_mem->hmin / SUNRabs(hcur));

  /* enforce maximum time step size */
  ark_mem->eta /= SUNMAX(ONE, SUNRabs(hcur) * ark_mem->hmax_inv * ark_mem->eta);

  SUNLogDebug(ARK_LOGGER, "new-step-eta", "eta = " SUN_FORMAT_G, ark_mem->eta);

  return (retval);
}

/*===============================================================
  EOF
  ===============================================================*/
