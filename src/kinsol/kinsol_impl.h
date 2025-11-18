/*
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * KINSOL solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _KINSOL_IMPL_H
#define _KINSOL_IMPL_H

#include <stdarg.h>

#include <kinsol/kinsol.h>
#include <sundials/priv/sundials_context_impl.h>

#include "sundials_iterative_impl.h"
#include "sundials_logger_impl.h"
#include "sundials_macros.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *   M A I N    S O L V E R    M E M O R Y    B L O C K
 * =================================================================
 */

/* KINSOL default constants */

#define MXITER_DEFAULT     200
#define MXNBCF_DEFAULT     10
#define MSBSET_DEFAULT     10
#define MSBSET_SUB_DEFAULT 5

#define OMEGA_MIN SUN_RCONST(0.00001)
#define OMEGA_MAX SUN_RCONST(0.9)

/*=================================================================*/
/* Shortcuts                                                       */
/*=================================================================*/

#define KIN_PROFILER kin_mem->kin_sunctx->profiler
#define KIN_LOGGER   kin_mem->kin_sunctx->logger

/*
 * -----------------------------------------------------------------
 * Types : struct KINMemRec and struct *KINMem
 * -----------------------------------------------------------------
 * A variable declaration of type struct *KINMem denotes a
 * pointer to a data structure of type struct KINMemRec. The
 * KINMemRec structure contains numerous fields that must be
 * accessible by KINSOL solver module routines.
 * -----------------------------------------------------------------
 */

typedef struct KINMemRec
{
  SUNContext kin_sunctx;

  sunrealtype kin_uround; /* machine epsilon (or unit roundoff error)
                                 (defined in sundials_types.h)                */

  /* problem specification data */

  KINSysFn kin_func;         /* nonlinear system function implementation     */
  void* kin_user_data;       /* work space available to func routine         */
  sunrealtype kin_fnormtol;  /* stopping tolerance on L2-norm of function
                                  value                                        */
  sunrealtype kin_scsteptol; /* scaled step length tolerance                 */
  int kin_globalstrategy;    /* choices are KIN_NONE, KIN_LINESEARCH
                                  KIN_PICARD and KIN_FP                        */
  long int kin_mxiter;       /* maximum number of nonlinear iterations       */
  long int kin_msbset;       /* maximum number of nonlinear iterations that
                                  may be performed between calls to the
                                  linear solver setup routine (lsetup)         */
  long int kin_msbset_sub;   /* subinterval length for residual monitoring   */
  long int kin_mxnbcf;       /* maximum number of beta condition failures    */
  int kin_etaflag;           /* choices are KIN_ETACONSTANT, KIN_ETACHOICE1
                                  and KIN_ETACHOICE2                           */
  sunbooleantype kin_noMinEps; /* flag controlling whether or not the value
                                  of eps is bounded below                      */
  sunbooleantype kin_constraintsSet; /* flag indicating if constraints are being
                                     used                                      */
  sunbooleantype kin_jacCurrent;     /* flag indicating if the Jacobian info.
                                     used by the linear solver is current      */
  sunbooleantype kin_callForcingTerm; /* flag set if using either KIN_ETACHOICE1
                                      or KIN_ETACHOICE2                        */
  sunbooleantype kin_noResMon;        /* flag indicating if the nonlinear
                                       residual monitoring scheme should be
                                       used                                    */
  sunbooleantype kin_retry_nni;       /* flag indicating if nonlinear iteration
                                       should be retried (set by residual
                                       monitoring algorithm)                   */
  sunbooleantype kin_update_fnorm_sub; /* flag indicating if the fnorm associated
                                       with the subinterval needs to be
                                       updated (set by residual monitoring
                                       algorithm)                              */

  sunrealtype kin_mxnewtstep; /* maximum allowable scaled step length         */
  sunrealtype kin_mxnstepin;  /* input (or preset) value for mxnewtstep       */
  sunrealtype kin_sqrt_relfunc; /* relative error bound for func(u)             */
  sunrealtype kin_stepl;     /* scaled length of current step                */
  sunrealtype kin_stepmul;   /* step scaling factor                          */
  sunrealtype kin_eps;       /* current value of eps                         */
  sunrealtype kin_eta;       /* current value of eta                         */
  sunrealtype kin_eta_gamma; /* gamma value used in eta calculation
                                  (choice #2)                                  */
  sunrealtype kin_eta_alpha; /* alpha value used in eta calculation
                                  (choice #2)                                  */
  sunbooleantype kin_noInitSetup; /* flag controlling whether or not the KINSol
                                  routine makes an initial call to the
                                  linear solver setup routine (lsetup)         */
  sunrealtype kin_sthrsh;         /* threshold value for calling the linear
                                  solver setup routine                         */

  /* counters */

  long int kin_nni;         /* number of nonlinear iterations               */
  long int kin_nfe;         /* number of calls made to func routine         */
  long int kin_nnilset;     /* value of nni counter when the linear solver
                                  setup was last called                        */
  long int kin_nnilset_sub; /* value of nni counter when the linear solver
                                  setup was last called (subinterval)          */
  long int kin_nbcf;        /* number of times the beta-condition could not
                                  be met in KINLineSearch                      */
  long int kin_nbktrk;      /* number of backtracks performed by
                                  KINLineSearch                                */
  long int kin_ncscmx;      /* number of consecutive steps of size
                                  mxnewtstep taken                             */

  /* vectors */

  N_Vector kin_uu;     /* solution vector/current iterate (initially
                               contains initial guess, but holds approximate
                               solution upon completion if no errors occurred) */
  N_Vector kin_unew;   /* next iterate (unew = uu+pp)                     */
  N_Vector kin_fval;   /* vector containing result of nonlinear system
                               function evaluated at a given iterate
                               (fval = func(uu))                               */
  N_Vector kin_gval;   /* vector containing result of the fixed point
                               function evaluated at a given iterate;
                               used in KIN_PICARD strategy only.
                               (gval = uu - L^{-1}fval(uu))                    */
  N_Vector kin_uscale; /* iterate scaling vector                          */
  N_Vector kin_fscale; /* fval scaling vector                             */
  N_Vector kin_pp;     /* incremental change vector (pp = unew-uu)        */
  N_Vector kin_constraints; /* constraints vector                              */
  N_Vector kin_vtemp1; /* scratch vector #1                               */
  N_Vector kin_vtemp2; /* scratch vector #2                               */
  N_Vector kin_vtemp3; /* scratch vector #3                               */

  /* fixed point and Picard options */
  sunbooleantype kin_ret_newest; /* return the newest FP iteration     */
  sunbooleantype kin_damping;    /* flag to apply damping in FP/Picard */
  sunrealtype kin_beta;          /* damping parameter for FP/Picard    */

  /* space requirements for AA, Broyden and NLEN */
  N_Vector kin_fold_aa;    /* vector needed for AA, Broyden, and NLEN         */
  N_Vector kin_gold_aa;    /* vector needed for AA, Broyden, and NLEN         */
  N_Vector* kin_df_aa;     /* vector array needed for AA, Broyden, and NLEN   */
  N_Vector* kin_dg_aa;     /* vector array needed for AA, Broyden and NLEN    */
  N_Vector* kin_q_aa;      /* vector array needed for AA                      */
  sunrealtype kin_beta_aa; /* beta damping parameter for AA                   */
  sunrealtype* kin_gamma_aa; /* array of size maa used in AA                    */
  sunrealtype* kin_R_aa;   /* array of size maa*maa used in AA                */
  sunrealtype* kin_T_aa;   /* array of size maa*maa used in AA with ICWY MGS  */
  long int kin_m_aa;       /* parameter for AA, Broyden or NLEN               */
  long int kin_m_aa_alloc; /* depth (m) used for AA memory allocations */
  long int kin_delay_aa;   /* number of iterations to delay AA */
  long int kin_current_depth;  /* current Anderson acceleration space size */
  KINDampingFn kin_damping_fn; /* function to determine the damping factor */
  KINDepthFn kin_depth_fn;     /* function to determine the depth with AA */
  int kin_orth_aa;             /* parameter for AA determining orthogonalization
                                 routine
                                 0 - Modified Gram Schmidt (standard)
                                 1 - ICWY Modified Gram Schmidt (Bjorck)
                                 2 - CGS2 (Hernandez)
                                 3 - Delayed CGS2 (Hernandez)                    */
  long int kin_orth_aa_alloc; /* depth (m) used for orthogonalization memory allocations */
  SUNQRAddFn kin_qr_func; /* QRAdd function for AA orthogonalization         */
  SUNQRData kin_qr_data;  /* Additional parameters required for QRAdd routine
                                 set for AA                                      */
  sunbooleantype kin_damping_aa; /* flag to apply damping in AA                     */
  sunbooleantype kin_dot_prod_sb; /* use single buffer dot product */
  sunrealtype* kin_cv; /* scalar array for fused vector operations        */
  N_Vector* kin_Xv;    /* vector array for fused vector operations        */

  /* space requirements for vector storage */

  sunindextype kin_lrw1; /* number of sunrealtype-sized memory blocks needed
                               for a single N_Vector                           */
  sunindextype kin_liw1; /* number of int-sized memory blocks needed for
                               a single N_Vecotr                               */
  long int kin_lrw;      /* total number of sunrealtype-sized memory blocks
                               needed for all KINSOL work vectors              */
  long int kin_liw;      /* total number of int-sized memory blocks needed
                               for all KINSOL work vectors                     */

  /* linear solver data */

  /* function prototypes (pointers) */

  int (*kin_linit)(struct KINMemRec* kin_mem);

  int (*kin_lsetup)(struct KINMemRec* kin_mem);

  int (*kin_lsolve)(struct KINMemRec* kin_mem, N_Vector xx, N_Vector bb,
                    sunrealtype* sJpnorm, sunrealtype* sFdotJp);

  int (*kin_lfree)(struct KINMemRec* kin_mem);

  sunbooleantype kin_inexact_ls; /* flag set by the linear solver module
                                 (in linit) indicating whether this is an
                                 iterative linear solver (SUNTRUE), or a direct
                                 linear solver (SUNFALSE)                       */

  void* kin_lmem; /* pointer to linear solver memory block             */

  sunrealtype kin_fnorm; /* value of L2-norm of fscale*fval                   */
  sunrealtype kin_f1norm; /* f1norm = 0.5*(fnorm)^2                            */
  sunrealtype kin_sFdotJp; /* value of scaled F(u) vector (fscale*fval)
                             dotted with scaled J(u)*pp vector (set by lsolve) */
  sunrealtype kin_sJpnorm; /* value of L2-norm of fscale*(J(u)*pp)
                             (set by lsolve)                                   */

  sunrealtype kin_fnorm_sub; /* value of L2-norm of fscale*fval (subinterval)     */
  sunbooleantype kin_eval_omega; /* flag indicating that omega must be evaluated. */
  sunrealtype kin_omega; /* constant value for real scalar used in test to
                             determine if reduction of norm of nonlinear
                             residual is sufficient. Unless a valid constant
                             value is specified by the user, omega is estimated
                             from omega_min and omega_max at each iteration.    */
  sunrealtype kin_omega_min; /* lower bound on omega                               */
  sunrealtype kin_omega_max; /* upper bound on omega                               */

  /*
   * -----------------------------------------------------------------
   * Note: The KINLineSearch subroutine scales the values of the
   * variables sFdotJp and sJpnorm by a factor rl (lambda) that is
   * chosen by the line search algorithm such that the sclaed Newton
   * step satisfies the following conditions:
   *
   *  F(u_k+1) <= F(u_k) + alpha*(F(u_k)^T * J(u_k))*p*rl
   *
   *  F(u_k+1) >= F(u_k) + beta*(F(u_k)^T * J(u_k))*p*rl
   *
   * where alpha = 1.0e-4, beta = 0.9, u_k+1 = u_k + rl*p,
   * 0 < rl <= 1, J denotes the system Jacobian, and F represents
   * the nonliner system function.
   * -----------------------------------------------------------------
   */

  sunbooleantype kin_MallocDone; /* flag indicating if KINMalloc has been
                                 called yet                                    */

}* KINMem;

/*
 * =================================================================
 *   I N T E R F A C E   T O   L I N E A R   S O L V E R
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_linit)(KINMem kin_mem)
 * -----------------------------------------------------------------
 * kin_linit initializes solver-specific data structures (including
 * variables used as counters or for storing statistical information),
 * but system memory allocation should be done by the subroutine
 * that actually initializes the environment for liner solver
 * package. If the linear system is to be preconditioned, then the
 * variable setupNonNull (type sunbooleantype) should be set to SUNTRUE
 * (predefined constant) and the kin_lsetup routine should be
 * appropriately defined.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 * If the necessary variables have been successfully initialized,
 * then the kin_linit function should return 0 (zero). Otherwise,
 * the subroutine should indicate a failure has occurred by
 * returning a non-zero integer value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_lsetup)(KINMem kin_mem)
 * -----------------------------------------------------------------
 * kin_lsetup interfaces with the user-supplied pset subroutine (the
 * preconditioner setup routine), and updates relevant variable
 * values (see KINSpgmrSetup/KINSpbcgSetup). Simply stated, the
 * kin_lsetup routine prepares the linear solver for a subsequent
 * call to the user-supplied kin_lsolve function.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 * If successful, the kin_lsetup routine should return 0 (zero).
 * Otherwise it should return a non-zero value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_lsolve)(KINMem kin_mem, N_Vector xx,
 *                N_Vector bb, sunrealtype *sJpnorm, sunrealtype *sFdotJp)
 * -----------------------------------------------------------------
 * kin_lsolve interfaces with the subroutine implementing the
 * numerical method to be used to solve the linear system J*xx = bb,
 * and must increment the relevant counter variable values in
 * addition to computing certain values used by the global strategy
 * and forcing term routines (see KINInexactNewton, KINLineSearch,
 * KINForcingTerm, and KINSpgmrSolve/KINSpbcgSolve).
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 *  xx  vector (type N_Vector) set to initial guess by kin_lsolve
 *      routine prior to calling the linear solver, but which upon
 *      return contains an approximate solution of the linear
 *      system J*xx = bb, where J denotes the system Jacobian
 *
 *  bb  vector (type N_Vector) set to -func(u) (negative of the
 *      value of the system function evaluated at the current
 *      iterate) by KINLinSolDrv before kin_lsolve is called
 *
 *  sJpnorm  holds the value of the L2-norm (Euclidean norm) of
 *           fscale*(J(u)*pp) upon return
 *
 *  sFdotJp  holds the value of the scaled F(u) (fscale*F) dotted
 *           with the scaled J(u)*pp vector upon return
 *
 * If successful, the kin_lsolve routine should return 0 (zero).
 * Otherwise it should return a positive value if a re-evaluation
 * of the lsetup function could recover, or a negative value if
 * no such recovery is possible.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_lfree)(KINMem kin_mem)
 * -----------------------------------------------------------------
 * kin_lfree is called by KINFree and should free (deallocate) all
 * system memory resources allocated for the linear solver module
 * (see KINSpgmrFree/KINSpbcgFree).  It should return 0 upon
 * success, nonzero on failure.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 *   K I N S O L    I N T E R N A L   F U N C T I O N S
 * =================================================================
 */

/* High level error handler */

void KINProcessError(KINMem kin_mem, int error_code, int line, const char* func,
                     const char* file, const char* msgfmt, ...);

/* High level info handler */

void KINPrintInfo(KINMem kin_mem, int info_code, const char* module,
                  const char* fname, const char* msgfmt, ...);

/* Prototype of internal infoHandler function */

void KINInfoHandler(const char* module, const char* function, char* msg,
                    void* user_data);

/* Anderson acceleration utilities */
int KINInitAA(KINMem kin_mem);
void KINFreeAA(KINMem kin_mem);

/* Orthogonalization utilities */
int KINInitOrth(KINMem kin_mem);
void KINFreeOrth(KINMem kin_mem);

/*
 * =================================================================
 *   K I N S O L    E R R O R    M E S S A G E S
 * =================================================================
 */

#define MSG_MEM_FAIL    "A memory request failed."
#define MSG_NO_MEM      "kinsol_mem = NULL illegal."
#define MSG_NULL_SUNCTX "sunctx = NULL illegal."
#define MSG_BAD_NVECTOR "A required vector operation is not implemented."
#define MSG_FUNC_NULL   "func = NULL illegal."
#define MSG_NO_MALLOC   "Attempt to call before KINMalloc illegal."

#define MSG_BAD_MXITER      "Illegal value for mxiter."
#define MSG_BAD_MSBSET      "Illegal msbset < 0."
#define MSG_BAD_MSBSETSUB   "Illegal msbsetsub < 0."
#define MSG_BAD_ETACHOICE   "Illegal value for etachoice."
#define MSG_BAD_ETACONST    "eta out of range."
#define MSG_BAD_GAMMA       "gamma out of range."
#define MSG_BAD_ALPHA       "alpha out of range."
#define MSG_BAD_MXNEWTSTEP  "Illegal mxnewtstep < 0."
#define MSG_BAD_RELFUNC     "relfunc < 0 illegal."
#define MSG_BAD_FNORMTOL    "fnormtol < 0 illegal."
#define MSG_BAD_SCSTEPTOL   "scsteptol < 0 illegal."
#define MSG_BAD_MXNBCF      "mxbcf < 0 illegal."
#define MSG_BAD_CONSTRAINTS "Illegal values in constraints vector."
#define MSG_BAD_OMEGA       "scalars < 0 illegal."
#define MSG_BAD_MAA         "maa < 0 illegal."
#define MSG_BAD_ORTHAA      "Illegal value for orthaa."
#define MSG_ZERO_MAA        "maa = 0 illegal."

#define MSG_LSOLV_NO_MEM       "The linear solver memory pointer is NULL."
#define MSG_UU_NULL            "uu = NULL illegal."
#define MSG_BAD_GLSTRAT        "Illegal value for global strategy."
#define MSG_BAD_USCALE         "uscale = NULL illegal."
#define MSG_USCALE_NONPOSITIVE "uscale has nonpositive elements."
#define MSG_BAD_FSCALE         "fscale = NULL illegal."
#define MSG_FSCALE_NONPOSITIVE "fscale has nonpositive elements."
#define MSG_CONSTRAINTS_NOTOK \
  "Constraints not allowed with fixed point or Picard iterations"
#define MSG_INITIAL_CNSTRNT "Initial guess does NOT meet constraints."
#define MSG_LINIT_FAIL      "The linear solver's init routine failed."

#define MSG_SYSFUNC_FAILED \
  "The system function failed in an unrecoverable manner."
#define MSG_SYSFUNC_FIRST "The system function failed at the first call."
#define MSG_LSETUP_FAILED \
  "The linear solver's setup function failed in an unrecoverable manner."
#define MSG_LSOLVE_FAILED \
  "The linear solver's solve function failed in an unrecoverable manner."
#define MSG_LINSOLV_NO_RECOVERY                                              \
  "The linear solver's solve function failed recoverably, but the Jacobian " \
  "data is already current."
#define MSG_LINESEARCH_NONCONV                                            \
  "The line search algorithm was unable to find an iterate sufficiently " \
  "distinct from the current iterate."
#define MSG_LINESEARCH_BCFAIL                                               \
  "The line search algorithm was unable to satisfy the beta-condition for " \
  "nbcfails iterations."
#define MSG_MAXITER_REACHED \
  "The maximum number of iterations was reached before convergence."
#define MSG_MXNEWT_5X_EXCEEDED                                                \
  "Five consecutive steps have been taken that satisfy a scaled step length " \
  "test."
#define MSG_SYSFUNC_REPTD \
  "Unable to correct repeated recoverable system function errors."
#define MSG_NOL_FAIL                                                  \
  "Unable to find user's Linear Jacobian, which is required for the " \
  "KIN_PICARD Strategy"

/*
 * =================================================================
 *   K I N S O L    I N F O    M E S S A G E S
 * =================================================================
 */

#define INFO_IVAR   "%s = %d"
#define INFO_LIVAR  "%s = %ld"
#define INFO_RETVAL "Return value: %d"
#define INFO_ADJ    "no. of lambda adjustments = %ld"

#define INFO_RVAR   "%s = " SUN_FORMAT_G
#define INFO_NNI    "nni = %4ld, nfe = %6ld, fnorm = " SUN_FORMAT_G
#define INFO_TOL    "scsteptol = " SUN_FORMAT_G ", fnormtol = " SUN_FORMAT_G
#define INFO_FMAX   "scaled f norm (for stopping) = " SUN_FORMAT_G
#define INFO_PNORM  "pnorm = " SUN_FORMAT_E
#define INFO_PNORM1 "(ivio=1) pnorm = " SUN_FORMAT_E
#define INFO_FNORM  "fnorm(L2) = " SUN_FORMAT_E
#define INFO_LAM                                                  \
  "min_lam = " SUN_FORMAT_E ", f1norm = " SUN_FORMAT_E ", pnorm " \
  "= " SUN_FORMAT_E
#define INFO_ALPHA                                   \
  "fnorm = " SUN_FORMAT_E ", f1norm = " SUN_FORMAT_E \
  ", alpha_cond = " SUN_FORMAT_E ", lam = " SUN_FORMAT_E
#define INFO_BETA                                                 \
  "f1norm = " SUN_FORMAT_E ", beta_cond = " SUN_FORMAT_E ", lam " \
  "= " SUN_FORMAT_E
#define INFO_ALPHABETA                                    \
  "f1norm = " SUN_FORMAT_E ", alpha_cond = " SUN_FORMAT_E \
  ", beta_cond = " SUN_FORMAT_E ", lam = " SUN_FORMAT_E

#ifdef __cplusplus
}
#endif

#endif
