/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the implementation file for ARKODE's root-finding (in
 * time) utility.
 *--------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_impl.h"

/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeRootInit:

  ARKodeRootInit initializes a rootfinding problem to be solved
  during the integration of the ODE system.  It loads the root
  function pointer and the number of root functions, notifies
  ARKODE that the "fullrhs" function is required, and allocates
  workspace memory.  The return value is ARK_SUCCESS = 0 if no
  errors occurred, or a negative value otherwise.
  ---------------------------------------------------------------*/
int ARKodeRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  int i, nrt;

  /* unpack ark_mem */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  nrt     = (nrtfn < 0) ? 0 : nrtfn;

  /* Ensure that stepper provides fullrhs function */
  if (nrt > 0)
  {
    if (!(ark_mem->step_fullrhs))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_MISSING_FULLRHS);
      return ARK_ILL_INPUT;
    }

    if (!arkAllocVec(ark_mem, ark_mem->yn, &ark_mem->fn))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
  }

  /* If unallocated, allocate rootfinding structure, set defaults, update space */
  if (ark_mem->root_mem == NULL)
  {
    ark_mem->root_mem = (ARKodeRootMem)malloc(sizeof(struct ARKodeRootMemRec));
    if (ark_mem->root_mem == NULL)
    {
      arkProcessError(ark_mem, 0, __LINE__, __func__, __FILE__,
                      MSG_ARK_ARKMEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->root_mem->glo       = NULL;
    ark_mem->root_mem->ghi       = NULL;
    ark_mem->root_mem->grout     = NULL;
    ark_mem->root_mem->iroots    = NULL;
    ark_mem->root_mem->rootdir   = NULL;
    ark_mem->root_mem->gfun      = NULL;
    ark_mem->root_mem->nrtfn     = 0;
    ark_mem->root_mem->irfnd     = 0;
    ark_mem->root_mem->gactive   = NULL;
    ark_mem->root_mem->mxgnull   = 1;
    ark_mem->root_mem->root_data = ark_mem->user_data;

    ark_mem->lrw += ARK_ROOT_LRW;
    ark_mem->liw += ARK_ROOT_LIW;
  }

  /* If rerunning ARKodeRootInit() with a different number of root
     functions (changing number of gfun components), then free
     currently held memory resources */
  if ((nrt != ark_mem->root_mem->nrtfn) && (ark_mem->root_mem->nrtfn > 0))
  {
    free(ark_mem->root_mem->glo);
    ark_mem->root_mem->glo = NULL;
    free(ark_mem->root_mem->ghi);
    ark_mem->root_mem->ghi = NULL;
    free(ark_mem->root_mem->grout);
    ark_mem->root_mem->grout = NULL;
    free(ark_mem->root_mem->iroots);
    ark_mem->root_mem->iroots = NULL;
    free(ark_mem->root_mem->rootdir);
    ark_mem->root_mem->rootdir = NULL;
    free(ark_mem->root_mem->gactive);
    ark_mem->root_mem->gactive = NULL;

    ark_mem->lrw -= 3 * (ark_mem->root_mem->nrtfn);
    ark_mem->liw -= 3 * (ark_mem->root_mem->nrtfn);
  }

  /* If ARKodeRootInit() was called with nrtfn == 0, then set
     nrtfn to zero and gfun to NULL before returning */
  if (nrt == 0)
  {
    ark_mem->root_mem->nrtfn = nrt;
    ark_mem->root_mem->gfun  = NULL;
    return (ARK_SUCCESS);
  }

  /* If rerunning ARKodeRootInit() with the same number of root
     functions (not changing number of gfun components), then
     check if the root function argument has changed */
  /* If g != NULL then return as currently reserved memory
     resources will suffice */
  if (nrt == ark_mem->root_mem->nrtfn)
  {
    if (g != ark_mem->root_mem->gfun)
    {
      if (g == NULL)
      {
        free(ark_mem->root_mem->glo);
        ark_mem->root_mem->glo = NULL;
        free(ark_mem->root_mem->ghi);
        ark_mem->root_mem->ghi = NULL;
        free(ark_mem->root_mem->grout);
        ark_mem->root_mem->grout = NULL;
        free(ark_mem->root_mem->iroots);
        ark_mem->root_mem->iroots = NULL;
        free(ark_mem->root_mem->rootdir);
        ark_mem->root_mem->rootdir = NULL;
        free(ark_mem->root_mem->gactive);
        ark_mem->root_mem->gactive = NULL;

        ark_mem->lrw -= 3 * nrt;
        ark_mem->liw -= 3 * nrt;

        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        MSG_ARK_NULL_G);
        return (ARK_ILL_INPUT);
      }
      else
      {
        ark_mem->root_mem->gfun = g;
        return (ARK_SUCCESS);
      }
    }
    else { return (ARK_SUCCESS); }
  }

  /* Set variable values in ARKODE memory block */
  ark_mem->root_mem->nrtfn = nrt;
  if (g == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_G);
    return (ARK_ILL_INPUT);
  }
  else { ark_mem->root_mem->gfun = g; }

  /* Allocate necessary memory and return */
  ark_mem->root_mem->glo = NULL;
  ark_mem->root_mem->glo = (sunrealtype*)malloc(nrt * sizeof(sunrealtype));
  if (ark_mem->root_mem->glo == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MEM_FAIL);
    return (ARK_MEM_FAIL);
  }
  ark_mem->root_mem->ghi = NULL;
  ark_mem->root_mem->ghi = (sunrealtype*)malloc(nrt * sizeof(sunrealtype));
  if (ark_mem->root_mem->ghi == NULL)
  {
    free(ark_mem->root_mem->glo);
    ark_mem->root_mem->glo = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MEM_FAIL);
    return (ARK_MEM_FAIL);
  }
  ark_mem->root_mem->grout = NULL;
  ark_mem->root_mem->grout = (sunrealtype*)malloc(nrt * sizeof(sunrealtype));
  if (ark_mem->root_mem->grout == NULL)
  {
    free(ark_mem->root_mem->glo);
    ark_mem->root_mem->glo = NULL;
    free(ark_mem->root_mem->ghi);
    ark_mem->root_mem->ghi = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MEM_FAIL);
    return (ARK_MEM_FAIL);
  }
  ark_mem->root_mem->iroots = NULL;
  ark_mem->root_mem->iroots = (int*)malloc(nrt * sizeof(int));
  if (ark_mem->root_mem->iroots == NULL)
  {
    free(ark_mem->root_mem->glo);
    ark_mem->root_mem->glo = NULL;
    free(ark_mem->root_mem->ghi);
    ark_mem->root_mem->ghi = NULL;
    free(ark_mem->root_mem->grout);
    ark_mem->root_mem->grout = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MEM_FAIL);
    return (ARK_MEM_FAIL);
  }
  ark_mem->root_mem->rootdir = NULL;
  ark_mem->root_mem->rootdir = (int*)malloc(nrt * sizeof(int));
  if (ark_mem->root_mem->rootdir == NULL)
  {
    free(ark_mem->root_mem->glo);
    ark_mem->root_mem->glo = NULL;
    free(ark_mem->root_mem->ghi);
    ark_mem->root_mem->ghi = NULL;
    free(ark_mem->root_mem->grout);
    ark_mem->root_mem->grout = NULL;
    free(ark_mem->root_mem->iroots);
    ark_mem->root_mem->iroots = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MEM_FAIL);
    return (ARK_MEM_FAIL);
  }
  ark_mem->root_mem->gactive = NULL;
  ark_mem->root_mem->gactive =
    (sunbooleantype*)malloc(nrt * sizeof(sunbooleantype));
  if (ark_mem->root_mem->gactive == NULL)
  {
    free(ark_mem->root_mem->glo);
    ark_mem->root_mem->glo = NULL;
    free(ark_mem->root_mem->ghi);
    ark_mem->root_mem->ghi = NULL;
    free(ark_mem->root_mem->grout);
    ark_mem->root_mem->grout = NULL;
    free(ark_mem->root_mem->iroots);
    ark_mem->root_mem->iroots = NULL;
    free(ark_mem->root_mem->rootdir);
    ark_mem->root_mem->rootdir = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MEM_FAIL);
    return (ARK_MEM_FAIL);
  }

  /* Set default values for rootdir (both directions) */
  for (i = 0; i < nrt; i++) { ark_mem->root_mem->rootdir[i] = 0; }

  /* Set default values for gactive (all active) */
  for (i = 0; i < nrt; i++) { ark_mem->root_mem->gactive[i] = SUNTRUE; }

  ark_mem->lrw += 3 * nrt;
  ark_mem->liw += 3 * nrt;

  return (ARK_SUCCESS);
}

/*===============================================================
  Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkRootFree

  This routine frees all memory associated with ARKODE's
  rootfinding module.
  ---------------------------------------------------------------*/
int arkRootFree(void* arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->root_mem != NULL)
  {
    if (ark_mem->root_mem->nrtfn > 0)
    {
      free(ark_mem->root_mem->glo);
      ark_mem->root_mem->glo = NULL;
      free(ark_mem->root_mem->ghi);
      ark_mem->root_mem->ghi = NULL;
      free(ark_mem->root_mem->grout);
      ark_mem->root_mem->grout = NULL;
      free(ark_mem->root_mem->iroots);
      ark_mem->root_mem->iroots = NULL;
      free(ark_mem->root_mem->rootdir);
      ark_mem->root_mem->rootdir = NULL;
      free(ark_mem->root_mem->gactive);
      ark_mem->root_mem->gactive = NULL;
      ark_mem->lrw -= 3 * ark_mem->root_mem->nrtfn;
      ark_mem->liw -= 3 * ark_mem->root_mem->nrtfn;
    }
    free(ark_mem->root_mem);
    ark_mem->lrw -= ARK_ROOT_LRW;
    ark_mem->liw -= ARK_ROOT_LIW;
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkPrintRootMem

  This routine outputs the root-finding memory structure to a
  specified file pointer.
  ---------------------------------------------------------------*/
int arkPrintRootMem(void* arkode_mem, FILE* outfile)
{
  int i;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->root_mem != NULL)
  {
    fprintf(outfile, "ark_nrtfn = %i\n", ark_mem->root_mem->nrtfn);
    fprintf(outfile, "ark_nge = %li\n", ark_mem->root_mem->nge);
    if (ark_mem->root_mem->iroots != NULL)
    {
      for (i = 0; i < ark_mem->root_mem->nrtfn; i++)
      {
        fprintf(outfile, "ark_iroots[%i] = %i\n", i,
                ark_mem->root_mem->iroots[i]);
      }
    }
    if (ark_mem->root_mem->rootdir != NULL)
    {
      for (i = 0; i < ark_mem->root_mem->nrtfn; i++)
      {
        fprintf(outfile, "ark_rootdir[%i] = %i\n", i,
                ark_mem->root_mem->rootdir[i]);
      }
    }
    fprintf(outfile, "ark_irfnd = %i\n", ark_mem->root_mem->irfnd);
    fprintf(outfile, "ark_mxgnull = %i\n", ark_mem->root_mem->mxgnull);
    if (ark_mem->root_mem->gactive != NULL)
    {
      for (i = 0; i < ark_mem->root_mem->nrtfn; i++)
      {
        fprintf(outfile, "ark_gactive[%i] = %i\n", i,
                ark_mem->root_mem->gactive[i]);
      }
    }
    fprintf(outfile, "ark_tlo = " SUN_FORMAT_G "\n", ark_mem->root_mem->tlo);
    fprintf(outfile, "ark_thi = " SUN_FORMAT_G "\n", ark_mem->root_mem->thi);
    fprintf(outfile, "ark_trout = " SUN_FORMAT_G "\n", ark_mem->root_mem->trout);
    if (ark_mem->root_mem->glo != NULL)
    {
      for (i = 0; i < ark_mem->root_mem->nrtfn; i++)
      {
        fprintf(outfile, "ark_glo[%i] = " SUN_FORMAT_G "\n", i,
                ark_mem->root_mem->glo[i]);
      }
    }
    if (ark_mem->root_mem->ghi != NULL)
    {
      for (i = 0; i < ark_mem->root_mem->nrtfn; i++)
      {
        fprintf(outfile, "ark_ghi[%i] = " SUN_FORMAT_G "\n", i,
                ark_mem->root_mem->ghi[i]);
      }
    }
    if (ark_mem->root_mem->grout != NULL)
    {
      for (i = 0; i < ark_mem->root_mem->nrtfn; i++)
      {
        fprintf(outfile, "ark_grout[%i] = " SUN_FORMAT_G "\n", i,
                ark_mem->root_mem->grout[i]);
      }
    }
    fprintf(outfile, "ark_ttol = " SUN_FORMAT_G "\n", ark_mem->root_mem->ttol);
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkRootCheck1

  This routine completes the initialization of rootfinding memory
  information, and checks whether g has a zero both at and very near
  the initial point of the IVP.

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0  if the g function failed, or
    ARK_SUCCESS     = 0  otherwise.
  ---------------------------------------------------------------*/
int arkRootCheck1(void* arkode_mem)
{
  int i, retval;
  sunrealtype smallh, hratio, tplus;
  sunbooleantype zroot;
  ARKodeMem ark_mem;
  ARKodeRootMem rootmem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  rootmem = ark_mem->root_mem;

  for (i = 0; i < rootmem->nrtfn; i++) { rootmem->iroots[i] = 0; }
  rootmem->tlo  = ark_mem->tcur;
  rootmem->ttol = (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h)) *
                  ark_mem->uround * HUND;

  /* Evaluate g at initial t and check for zero values. */
  retval       = rootmem->gfun(rootmem->tlo, ark_mem->yn, rootmem->glo,
                               rootmem->root_data);
  rootmem->nge = 1;
  if (retval != 0)
  {
    arkProcessError(ark_mem, ARK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
    return (ARK_RTFUNC_FAIL);
  }

  zroot = SUNFALSE;
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    if (SUNRabs(rootmem->glo[i]) == ZERO)
    {
      zroot               = SUNTRUE;
      rootmem->gactive[i] = SUNFALSE;
    }
  }
  if (!zroot) { return (ARK_SUCCESS); }

  /* call full RHS if needed */
  if (!(ark_mem->fn_is_current))
  {
    retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tn, ark_mem->yn,
                                   ark_mem->fn, ARK_FULLRHS_START);
    if (retval)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, ark_mem->tcur);
      return ARK_RHSFUNC_FAIL;
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  hratio = SUNMAX(rootmem->ttol / SUNRabs(ark_mem->h), TENTH);
  smallh = hratio * ark_mem->h;
  tplus  = rootmem->tlo + smallh;
  N_VLinearSum(ONE, ark_mem->yn, smallh, ark_mem->fn, ark_mem->ycur);
  retval = rootmem->gfun(tplus, ark_mem->ycur, rootmem->ghi, rootmem->root_data);
  rootmem->nge++;
  if (retval != 0)
  {
    arkProcessError(ark_mem, ARK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
    return (ARK_RTFUNC_FAIL);
  }

  /* We check now only the components of g which were exactly 0.0 at t0
   * to see if we can 'activate' them. */
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    if (!rootmem->gactive[i] && SUNRabs(rootmem->ghi[i]) != ZERO)
    {
      rootmem->gactive[i] = SUNTRUE;
      rootmem->glo[i]     = rootmem->ghi[i];
    }
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkRootCheck2

  This routine checks for exact zeros of g at the last root found,
  if the last return was a root.  It then checks for a close pair of
  zeros (an error condition), and for a new root at a nearby point.
  The array glo = g(tlo) at the left endpoint of the search interval
  is adjusted if necessary to assure that all g_i are nonzero
  there, before returning to do a root search in the interval.

  On entry, tlo = tretlast is the last value of tret returned by
  ARKODE.  This may be the previous tn, the previous tout value, or
  the last root location.

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0 if the g function failed, or
    CLOSERT         = 3 if a close pair of zeros was found, or
    RTFOUND         = 1 if a new zero of g was found near tlo, or
    ARK_SUCCESS     = 0 otherwise.
  ---------------------------------------------------------------*/
int arkRootCheck2(void* arkode_mem)
{
  int i, retval;
  sunrealtype smallh, tplus;
  sunbooleantype zroot;
  ARKodeMem ark_mem;
  ARKodeRootMem rootmem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  rootmem = ark_mem->root_mem;

  /* return if no roots in previous step */
  if (rootmem->irfnd == 0) { return (ARK_SUCCESS); }

  /* Set ark_ycur = y(tlo) */
  (void)ARKodeGetDky(ark_mem, rootmem->tlo, 0, ark_mem->ycur);

  /* Evaluate root-finding function: glo = g(tlo, y(tlo)) */
  retval = rootmem->gfun(rootmem->tlo, ark_mem->ycur, rootmem->glo,
                         rootmem->root_data);
  rootmem->nge++;
  if (retval != 0) { return (ARK_RTFUNC_FAIL); }

  /* reset root-finding flags (overall, and for specific eqns) */
  zroot = SUNFALSE;
  for (i = 0; i < rootmem->nrtfn; i++) { rootmem->iroots[i] = 0; }

  /* for all active roots, check if glo_i == 0 to mark roots found */
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    if (!rootmem->gactive[i]) { continue; }
    if (SUNRabs(rootmem->glo[i]) == ZERO)
    {
      zroot              = SUNTRUE;
      rootmem->iroots[i] = 1;
    }
  }
  if (!zroot) { return (ARK_SUCCESS); /* return if no roots */ }

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  /*     set time tolerance */
  rootmem->ttol = (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h)) *
                  ark_mem->uround * HUND;
  /*     set tplus = tlo + smallh */
  smallh = (ark_mem->h > ZERO) ? rootmem->ttol : -rootmem->ttol;
  tplus  = rootmem->tlo + smallh;
  /*     update ark_ycur with small explicit Euler step (if tplus is past tn) */
  if ((tplus - ark_mem->tcur) * ark_mem->h >= ZERO)
  {
    /* hratio = smallh/ark_mem->h; */
    N_VLinearSum(ONE, ark_mem->ycur, smallh, ark_mem->fn, ark_mem->ycur);
  }
  else
  {
    /*   set ark_ycur = y(tplus) via interpolation */
    (void)ARKodeGetDky(ark_mem, tplus, 0, ark_mem->ycur);
  }
  /*     set ghi = g(tplus,y(tplus)) */
  retval = rootmem->gfun(tplus, ark_mem->ycur, rootmem->ghi, rootmem->root_data);
  rootmem->nge++;
  if (retval != 0) { return (ARK_RTFUNC_FAIL); }

  /* Check for close roots (error return), for a new zero at tlo+smallh,
  and for a g_i that changed from zero to nonzero. */
  zroot = SUNFALSE;
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    if (!rootmem->gactive[i]) { continue; }
    if (SUNRabs(rootmem->ghi[i]) == ZERO)
    {
      if (rootmem->iroots[i] == 1) { return (CLOSERT); }
      zroot              = SUNTRUE;
      rootmem->iroots[i] = 1;
    }
    else
    {
      if (rootmem->iroots[i] == 1) { rootmem->glo[i] = rootmem->ghi[i]; }
    }
  }
  if (zroot) { return (RTFOUND); }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkRootCheck3

  This routine interfaces to arkRootfind to look for a root of g
  between tlo and either tn or tout, whichever comes first.
  Only roots beyond tlo in the direction of integration are sought.

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0 if the g function failed, or
    RTFOUND         = 1 if a root of g was found, or
    ARK_SUCCESS     = 0 otherwise.
  ---------------------------------------------------------------*/
int arkRootCheck3(void* arkode_mem, sunrealtype tout, int itask)
{
  int i, retval, ier;
  ARKodeMem ark_mem;
  ARKodeRootMem rootmem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  rootmem = ark_mem->root_mem;

  /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
  if (itask == ARK_ONE_STEP)
  {
    rootmem->thi = ark_mem->tcur;
    N_VScale(ONE, ark_mem->yn, ark_mem->ycur);
  }
  if (itask == ARK_NORMAL)
  {
    if ((tout - ark_mem->tcur) * ark_mem->h >= ZERO)
    {
      rootmem->thi = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, ark_mem->ycur);
    }
    else
    {
      rootmem->thi = tout;
      (void)ARKodeGetDky(ark_mem, rootmem->thi, 0, ark_mem->ycur);
    }
  }

  /* Set rootmem->ghi = g(thi) and call arkRootfind to search (tlo,thi) for roots. */
  retval = rootmem->gfun(rootmem->thi, ark_mem->ycur, rootmem->ghi,
                         rootmem->root_data);
  rootmem->nge++;
  if (retval != 0) { return (ARK_RTFUNC_FAIL); }

  rootmem->ttol = (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h)) *
                  ark_mem->uround * HUND;
  ier = arkRootfind(ark_mem);
  if (ier == ARK_RTFUNC_FAIL) { return (ARK_RTFUNC_FAIL); }
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    if (!rootmem->gactive[i] && rootmem->grout[i] != ZERO)
    {
      rootmem->gactive[i] = SUNTRUE;
    }
  }
  rootmem->tlo = rootmem->trout;
  for (i = 0; i < rootmem->nrtfn; i++) { rootmem->glo[i] = rootmem->grout[i]; }

  /* If no root found, return ARK_SUCCESS. */
  if (ier == ARK_SUCCESS) { return (ARK_SUCCESS); }

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void)ARKodeGetDky(ark_mem, rootmem->trout, 0, ark_mem->ycur);
  return (RTFOUND);
}

/*---------------------------------------------------------------
  arkRootfind

  This routine solves for a root of g(t) between tlo and thi, if
  one exists.  Only roots of odd multiplicity (i.e. with a change
  of sign in one of the g_i), or exact zeros, are found.
  Here the sign of tlo - thi is arbitrary, but if multiple roots
  are found, the one closest to tlo is returned.

  The method used is the Illinois algorithm, a modified secant method.
  Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
  Defined Output Points for Solutions of ODEs, Sandia National
  Laboratory Report SAND80-0180, February 1980.

  This routine uses the following parameters for communication:

  nrtfn    = number of functions g_i, or number of components of
            the vector-valued function g(t).  Input only.

  gfun     = user-defined function for g(t).  Its form is
             (void) gfun(t, y, gt, user_data)

  rootdir  = in array specifying the direction of zero-crossings.
             If rootdir[i] > 0, search for roots of g_i only if
             g_i is increasing; if rootdir[i] < 0, search for
             roots of g_i only if g_i is decreasing; otherwise
             always search for roots of g_i.

  gactive  = array specifying whether a component of g should
             or should not be monitored. gactive[i] is initially
             set to SUNTRUE for all i=0,...,nrtfn-1, but it may be
             reset to SUNFALSE if at the first step g[i] is 0.0
             both at the I.C. and at a small perturbation of them.
             gactive[i] is then set back on SUNTRUE only after the
             corresponding g function moves away from 0.0.

  nge      = cumulative counter for gfun calls.

  ttol     = a convergence tolerance for trout.  Input only.
             When a root at trout is found, it is located only to
             within a tolerance of ttol.  Typically, ttol should
             be set to a value on the order of
               100 * UROUND * max (SUNRabs(tlo), SUNRabs(thi))
             where UROUND is the unit roundoff of the machine.

  tlo, thi = endpoints of the interval in which roots are sought.
             On input, and must be distinct, but tlo - thi may
             be of either sign.  The direction of integration is
             assumed to be from tlo to thi.  On return, tlo and thi
             are the endpoints of the final relevant interval.

  glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
             and g(thi) respectively.  Input and output.  On input,
             none of the glo[i] should be zero.

  trout    = root location, if a root was found, or thi if not.
             Output only.  If a root was found other than an exact
             zero of g, trout is the endpoint thi of the final
             interval bracketing the root, with size at most ttol.

  grout    = array of length nrtfn containing g(trout) on return.

  iroots   = int array of length nrtfn with root information.
             Output only.  If a root was found, iroots indicates
             which components g_i have a root at trout.  For
             i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
             and g_i is increasing, iroots[i] = -1 if g_i has a
             root and g_i is decreasing, and iroots[i] = 0 if g_i
             has no roots or g_i varies in the direction opposite
             to that indicated by rootdir[i].

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0 if the g function failed, or
    RTFOUND         = 1 if a root of g was found, or
    ARK_SUCCESS     = 0 otherwise.
  ---------------------------------------------------------------*/
int arkRootfind(void* arkode_mem)
{
  sunrealtype alpha, tmid, gfrac, maxfrac, fracint, fracsub;
  int i, retval, imax, side, sideprev;
  sunbooleantype zroot, sgnchg;
  ARKodeMem ark_mem;
  ARKodeRootMem rootmem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  rootmem = ark_mem->root_mem;

  imax = 0;

  /* First check for change in sign in ghi or for a zero in ghi. */
  maxfrac = ZERO;
  zroot   = SUNFALSE;
  sgnchg  = SUNFALSE;
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    if (!rootmem->gactive[i]) { continue; }
    if (SUNRabs(rootmem->ghi[i]) == ZERO)
    {
      if (rootmem->rootdir[i] * rootmem->glo[i] <= ZERO) { zroot = SUNTRUE; }
    }
    else
    {
      if ((SUNRdifferentsign(rootmem->glo[i], rootmem->ghi[i])) &&
          (rootmem->rootdir[i] * rootmem->glo[i] <= ZERO))
      {
        gfrac = SUNRabs(rootmem->ghi[i] / (rootmem->ghi[i] - rootmem->glo[i]));
        if (gfrac > maxfrac)
        {
          sgnchg  = SUNTRUE;
          maxfrac = gfrac;
          imax    = i;
        }
      }
    }
  }

  /* If no sign change was found, reset trout and grout.  Then return
     ARK_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */
  if (!sgnchg)
  {
    rootmem->trout = rootmem->thi;
    for (i = 0; i < rootmem->nrtfn; i++)
    {
      rootmem->grout[i] = rootmem->ghi[i];
    }
    if (!zroot) { return (ARK_SUCCESS); }
    for (i = 0; i < rootmem->nrtfn; i++)
    {
      rootmem->iroots[i] = 0;
      if (!rootmem->gactive[i]) { continue; }
      if (SUNRabs(rootmem->ghi[i]) == ZERO)
      {
        rootmem->iroots[i] = rootmem->glo[i] > 0 ? -1 : 1;
      }
    }
    return (RTFOUND);
  }

  /* Initialize alpha to avoid compiler warning */
  alpha = ONE;

  /* A sign change was found.  Loop to locate nearest root. */
  side     = 0;
  sideprev = -1;
  for (;;)
  { /* Looping point */

    /* If interval size is already less than tolerance ttol, break. */
    if (SUNRabs(rootmem->thi - rootmem->tlo) <= rootmem->ttol) { break; }

    /* Set weight alpha.
       On the first two passes, set alpha = 1.  Thereafter, reset alpha
       according to the side (low vs high) of the subinterval in which
       the sign change was found in the previous two passes.
       If the sides were opposite, set alpha = 1.
       If the sides were the same, then double alpha (if high side),
       or halve alpha (if low side).
       The next guess tmid is the secant method value if alpha = 1, but
       is closer to tlo if alpha < 1, and closer to thi if alpha > 1.    */
    if (sideprev == side) { alpha = (side == 2) ? alpha * TWO : alpha * HALF; }
    else { alpha = ONE; }

    /* Set next root approximation tmid and get g(tmid).
       If tmid is too close to tlo or thi, adjust it inward,
       by a fractional distance that is between 0.1 and 0.5.  */
    tmid = rootmem->thi - (rootmem->thi - rootmem->tlo) * rootmem->ghi[imax] /
                            (rootmem->ghi[imax] - alpha * rootmem->glo[imax]);
    if (SUNRabs(tmid - rootmem->tlo) < HALF * rootmem->ttol)
    {
      fracint = SUNRabs(rootmem->thi - rootmem->tlo) / rootmem->ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF / fracint;
      tmid    = rootmem->tlo + fracsub * (rootmem->thi - rootmem->tlo);
    }
    if (SUNRabs(rootmem->thi - tmid) < HALF * rootmem->ttol)
    {
      fracint = SUNRabs(rootmem->thi - rootmem->tlo) / rootmem->ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF / fracint;
      tmid    = rootmem->thi - fracsub * (rootmem->thi - rootmem->tlo);
    }

    (void)ARKodeGetDky(ark_mem, tmid, 0, ark_mem->ycur);
    retval = rootmem->gfun(tmid, ark_mem->ycur, rootmem->grout,
                           rootmem->root_data);
    rootmem->nge++;
    if (retval != 0) { return (ARK_RTFUNC_FAIL); }

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */
    maxfrac  = ZERO;
    zroot    = SUNFALSE;
    sgnchg   = SUNFALSE;
    sideprev = side;
    for (i = 0; i < rootmem->nrtfn; i++)
    {
      if (!rootmem->gactive[i]) { continue; }
      if (SUNRabs(rootmem->grout[i]) == ZERO)
      {
        if (rootmem->rootdir[i] * rootmem->glo[i] <= ZERO) { zroot = SUNTRUE; }
      }
      else
      {
        if ((SUNRdifferentsign(rootmem->glo[i], rootmem->grout[i])) &&
            (rootmem->rootdir[i] * rootmem->glo[i] <= ZERO))
        {
          gfrac =
            SUNRabs(rootmem->grout[i] / (rootmem->grout[i] - rootmem->glo[i]));
          if (gfrac > maxfrac)
          {
            sgnchg  = SUNTRUE;
            maxfrac = gfrac;
            imax    = i;
          }
        }
      }
    }
    if (sgnchg)
    {
      /* Sign change found in (tlo,tmid); replace thi with tmid. */
      rootmem->thi = tmid;
      for (i = 0; i < rootmem->nrtfn; i++)
      {
        rootmem->ghi[i] = rootmem->grout[i];
      }
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (SUNRabs(rootmem->thi - rootmem->tlo) <= rootmem->ttol) { break; }
      continue; /* Return to looping point. */
    }

    if (zroot)
    {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      rootmem->thi = tmid;
      for (i = 0; i < rootmem->nrtfn; i++)
      {
        rootmem->ghi[i] = rootmem->grout[i];
      }
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    rootmem->tlo = tmid;
    for (i = 0; i < rootmem->nrtfn; i++)
    {
      rootmem->glo[i] = rootmem->grout[i];
    }
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (SUNRabs(rootmem->thi - rootmem->tlo) <= rootmem->ttol) { break; }

  } /* End of root-search loop */

  /* Reset trout and grout, set iroots, and return RTFOUND. */
  rootmem->trout = rootmem->thi;
  for (i = 0; i < rootmem->nrtfn; i++)
  {
    rootmem->grout[i]  = rootmem->ghi[i];
    rootmem->iroots[i] = 0;
    if (!rootmem->gactive[i]) { continue; }
    if ((SUNRabs(rootmem->ghi[i]) == ZERO) &&
        (rootmem->rootdir[i] * rootmem->glo[i] <= ZERO))
    {
      rootmem->iroots[i] = rootmem->glo[i] > 0 ? -1 : 1;
    }
    if ((SUNRdifferentsign(rootmem->glo[i], rootmem->ghi[i])) &&
        (rootmem->rootdir[i] * rootmem->glo[i] <= ZERO))
    {
      rootmem->iroots[i] = rootmem->glo[i] > 0 ? -1 : 1;
    }
  }
  return (RTFOUND);
}

/*===============================================================
  EOF
  ===============================================================*/
