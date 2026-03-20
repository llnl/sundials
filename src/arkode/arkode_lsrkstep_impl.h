/*---------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ UMBC
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
 * Implementation header file for ARKODE's LSRK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_LSRKSTEP_IMPL_H
#define _ARKODE_LSRKSTEP_IMPL_H

#include <arkode/arkode_lsrkstep.h>
#include <sundials/sundials_domeigestimator.h>

#include "arkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

#define STAGE_MAX_LIMIT_DEFAULT          200
#define DOM_EIG_SAFETY_DEFAULT           SUN_RCONST(1.01)
#define DOM_EIG_FREQ_DEFAULT             25
#define DOM_EIG_NUM_WARMUPS_DEFAULT      0
#define DOM_EIG_NUM_INIT_WARMUPS_DEFAULT -1 /* use DEE's default value */

/*===============================================================
  LSRK time step module private math function macros
  ===============================================================
 * SUNRlog calls the appropriate version of log
 *
 * SUNRsinh calls the appropriate version of sinh
 *
 * SUNRcosh calls the appropriate version of cosh
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : SUNRlog
 * -----------------------------------------------------------------
 * Usage : sunrealtype log_x;
 *         log_x = SUNRlog(x);
 * -----------------------------------------------------------------
 * SUNRlog(x) returns log(x) (base-e logarithmic function).
 * -----------------------------------------------------------------
 */

#ifndef SUNRlog
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRlog(x) (log((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRlog(x) (logf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRlog(x) (logl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRsinh
 * -----------------------------------------------------------------
 * Usage : sunrealtype sinh_x;
 *         sinh_x = SUNRsinh(x);
 * -----------------------------------------------------------------
 * SUNRsinh(x) returns sinh(x) (the hyperbolic sine of x).
 * -----------------------------------------------------------------
 */

#ifndef SUNRsinh
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsinh(x) (sinh((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsinh(x) (sinhf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsinh(x) (sinhl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*
 * -----------------------------------------------------------------
 * Function : SUNRcosh
 * -----------------------------------------------------------------
 * Usage : sunrealtype cosh_x;
 *         cosh_x = SUNRcosh(x);
 * -----------------------------------------------------------------
 * SUNRcosh(x) returns cosh(x) (the hyperbolic cosine of x).
 * -----------------------------------------------------------------
 */

#ifndef SUNRcosh
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRcosh(x) (cosh((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRcosh(x) (coshf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRcosh(x) (coshl((x)))
#else
#error \
  "SUNDIALS precision not defined, report to github.com/LLNL/sundials/issues"
#endif
#endif

/*===============================================================
  LSRK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeLSRKStepMemRec, ARKodeLSRKStepMem
  ---------------------------------------------------------------
  The type ARKodeLSRKStepMem is type pointer to struct
  ARKodeLSRKStepMemRec.  This structure contains fields to
  perform an explicit Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeLSRKStepMemRec
{
  /* LSRK problem specification */
  ARKRhsFn fe;
  ARKDomEigFn dom_eig_fn;

  int q; /* method order               */
  int p; /* embedding order            */

  int req_stages; /* number of requested stages   */

  ARKODE_LSRKMethodType LSRKmethod;

  /* Counters and stats*/
  long int nfe;   /* num fe calls       */
  long int nfeDQ; /* num fe calls for difference quotient approximation */
  long int dom_eig_num_evals; /* num of dom_eig computations   */
  int stage_max;              /* num of max stages used      */
  int stage_max_limit;        /* max allowed num of stages     */
  long int dom_eig_nst; /* num of step at which the last domainant eigenvalue was computed  */
  long int step_nst;      /* The number of successful steps. */
  long int num_dee_iters; /* number of iterations in the DEE estimates */

  /* Spectral info */
  sunrealtype lambdaR;         /* Real part of the dominated eigenvalue*/
  sunrealtype lambdaI;         /* Imaginary part of the dominated eigenvalue*/
  sunrealtype spectral_radius; /* spectral radius*/
  sunrealtype spectral_radius_max; /* max spectral radius*/
  sunrealtype spectral_radius_min; /* min spectral radius*/
  sunrealtype dom_eig_safety; /* some safety factor for the user provided dom_eig*/
  long int dom_eig_freq; /* indicates dom_eig update after dom_eig_freq successful steps*/
  int num_init_warmups; /* number of warm-ups in the first DEE estimates */
  int num_warmups;      /* number of warm-ups in succeeding DEE estimates */

  SUNDomEigEstimator DEE; /* DomEig estimator*/

  /* Flags */
  sunbooleantype dom_eig_update; /* flag indicating new dom_eig is needed */
  sunbooleantype const_Jac;      /* flag indicating Jacobian is constant */
  sunbooleantype dom_eig_is_current; /* SUNTRUE if dom_eig has been evaluated at tn */
  sunbooleantype is_SSP;             /* flag indicating SSP method*/
  sunbooleantype init_warmup;        /* flag indicating initial warm-up*/

  /* Reusable fused vector operation arrays */
  sunrealtype* cvals;
  N_Vector* Xvecs;
  int nfusedopvecs; /* length of cvals and Xvecs arrays */

}* ARKodeLSRKStepMem;

/*===============================================================
  LSRK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
void* lsrkStep_Create_Commons(ARKRhsFn rhs, sunrealtype t0, N_Vector y0,
                              SUNContext sunctx);
int lsrkStep_ReInit_Commons(void* arkode_mem, ARKRhsFn rhs, sunrealtype t0,
                            N_Vector y0);
int lsrkStep_Init(ARKodeMem ark_mem, sunrealtype tout, int init_type);
int lsrkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode);
int lsrkStep_TakeStepRKC(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepRKL(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSPs2(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSPs3(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSP43(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSP104(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                            int* nflagPtr);
int lsrkStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                        size_t offset, sunbooleantype* arg_used);
int lsrkStep_SetDefaults(ARKodeMem ark_mem);
int lsrkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int lsrkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
void lsrkStep_Free(ARKodeMem ark_mem);
void lsrkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);
int lsrkStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                            long int* rhs_evals);
int lsrkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele);

/* Internal utility routines */
int lsrkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem, ARKodeLSRKStepMem* step_mem);
int lsrkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodeLSRKStepMem* step_mem);
void lsrkStep_DomEigUpdateLogic(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem,
                                sunrealtype dsm);
int lsrkStep_ComputeNewDomEig(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem);
int lsrkStep_DQJtimes(void* arkode_mem, N_Vector v, N_Vector Jv);

/*===============================================================
  Reusable LSRKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_LSRKSTEP_NO_MEM "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
