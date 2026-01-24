/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 * Implementation header file for the main CVODES integrator.
 * -----------------------------------------------------------------*/

#ifndef _CVODES_IMPL_H
#define _CVODES_IMPL_H

#include <stdarg.h>

#include <cvodes/cvodes.h>
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/sundials_math.h>

#include "cvodes_proj_impl.h"
#include "sundials_logger_impl.h"
#include "sundials_macros.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================*/
/* Shortcuts                                                       */
/*=================================================================*/

#define CV_PROFILER cv_mem->cv_sunctx->profiler
#define CV_LOGGER   cv_mem->cv_sunctx->logger

/*
 * =================================================================
 *   I N T E R N A L   C O N S T A N T S
 * =================================================================
 */

/* Basic constants */

#define ADAMS_Q_MAX 12          /* max value of q for lmm == ADAMS     */
#define BDF_Q_MAX   5           /* max value of q for lmm == BDF       */
#define Q_MAX       ADAMS_Q_MAX /* max value of q for either lmm       */
#define L_MAX       (Q_MAX + 1) /* max value of L for either lmm       */
#define NUM_TESTS   5           /* number of error test quantities     */

#define HMIN_DEFAULT     SUN_RCONST(0.0) /* hmin default value     */
#define HMAX_INV_DEFAULT SUN_RCONST(0.0) /* hmax_inv default value */
#define MXHNIL_DEFAULT   10              /* mxhnil default value   */
#define MXSTEP_DEFAULT   500             /* mxstep default value   */

#define MSBP_DEFAULT 20 /* max steps between lsetup calls */
#define DGMAX_LSETUP_DEFAULT \
  SUN_RCONST(0.3) /* gamma threshold to call lsetup */

/* Step size change constants
 * --------------------------
 * ETA_MIN_FX_DEFAULT  if eta_min_fx < eta < eta_max_fx reject a change in step size or order
 * ETA_MAX_FX_DEFAULT
 * ETA_MAX_FS_DEFAULT  -+
 * ETA_MAX_ES_DEFAULT   |
 * ETA_MAX_GS_DEFAULT   |
 * ETA_MIN_DEFAULT      |-> bounds for eta (step size change)
 * ETA_MAX_EF_DEFAULT   |
 * ETA_MIN_EF_DEFAULT   |
 * ETA_CF_DEFAULT      -+
 * SMALL_NST_DEFAULT   nst <= SMALL_NST => use eta_max_es
 * SMALL_NEF_DEFAULT   if small_nef <= nef <= MXNEF1, then eta =  SUNMIN(eta, eta_max_ef)
 * ONEPSM (1+epsilon)  used in testing if the step size is below its bound
 */

#define ETA_MIN_FX_DEFAULT SUN_RCONST(0.0)
#define ETA_MAX_FX_DEFAULT SUN_RCONST(1.5)
#define ETA_MAX_FS_DEFAULT SUN_RCONST(10000.0)
#define ETA_MAX_ES_DEFAULT SUN_RCONST(10.0)
#define ETA_MAX_GS_DEFAULT SUN_RCONST(10.0)
#define ETA_MIN_DEFAULT    SUN_RCONST(0.1)
#define ETA_MAX_EF_DEFAULT SUN_RCONST(0.2)
#define ETA_MIN_EF_DEFAULT SUN_RCONST(0.1)
#define ETA_CF_DEFAULT     SUN_RCONST(0.25)
#define SMALL_NST_DEFAULT  10
#define SMALL_NEF_DEFAULT  2
#define ONEPSM             SUN_RCONST(1.000001)

/* Step size controller constants
 * ------------------------------
 * ADDON  safety factor in computing eta
 * BIAS1  -+
 * BIAS2   |-> bias factors in eta selection
 * BIAS3  -+
 */

#define ADDON SUN_RCONST(0.000001)
#define BIAS1 SUN_RCONST(6.0)
#define BIAS2 SUN_RCONST(6.0)
#define BIAS3 SUN_RCONST(10.0)

/* Order selection constants
 * -------------------------
 * LONG_WAIT   number of steps to wait before considering an order change when
 *             q==1 and MXNEF1 error test failures have occurred
 */

#define LONG_WAIT 10

/* Failure limits
 * --------------
 * MXNCF   max no. of convergence failures during one step try
 * MXNEF   max no. of error test failures during one step try
 * MXNEF1  max no. of error test failures before forcing a reduction of order
 */

#define MXNCF  10
#define MXNEF  7
#define MXNEF1 3

/* Control constants for lower-level functions used by cvStep
 * ----------------------------------------------------------
 *
 * cvHin return values:
 *    CV_SUCCESS,
 *    CV_RHSFUNC_FAIL,  CV_RPTD_RHSFUNC_ERR,
 *    CV_QRHSFUNC_FAIL, CV_RPTD_QRHSFUNC_ERR,
 *    CV_SRHSFUNC_FAIL, CV_RPTD_SRHSFUNC_ERR,
 *    CV_TOO_CLOSE
 *
 * cvStep control constants:
 *    DO_ERROR_TEST
 *    PREDICT_AGAIN
 *
 * cvStep return values:
 *    CV_SUCCESS,
 *    CV_CONV_FAILURE,      CV_ERR_FAILURE,
 *    CV_LSETUP_FAIL,       CV_LSOLVE_FAIL,
 *    CV_RTFUNC_FAIL,
 *    CV_RHSFUNC_FAIL,      CV_QRHSFUNC_FAIL,      CV_SRHSFUNC_FAIL,      CV_QSRHSFUNC_FAIL,
 *    CV_FIRST_RHSFUNC_ERR, CV_FIRST_QRHSFUNC_ERR, CV_FIRST_SRHSFUNC_ERR, CV_FIRST_QSRHSFUNC_ERR,
 *    CV_UNREC_RHSFUNC_ERR, CV_UNREC_QRHSFUNC_ERR, CV_UNREC_SRHSFUNC_ERR, CV_UNREC_QSRHSFUNC_ERR,
 *    CV_REPTD_RHSFUNC_ERR, CV_REPTD_QRHSFUNC_ERR, CV_REPTD_SRHSFUNC_ERR, CV_REPTD_QSRHSFUNC_ERR,
 *
 * cvNls input nflag values:
 *    FIRST_CALL
 *    PREV_CONV_FAIL
 *    PREV_PROJ_FAIL
 *    PREV_ERR_FAIL
 *
 * cvNls return values:
 *    CV_SUCCESS,
 *    CV_LSETUP_FAIL,     CV_LSOLVE_FAIL,
 *    CV_RHSFUNC_FAIL,    CV_SRHSFUNC_FAIL,
 *    SUN_NLS_CONV_RECVR,
 *    RHSFUNC_RECVR,      SRHSFUNC_RECVR
 *
 * cvNewtonIteration return values:
 *    CV_SUCCESS,
 *    CV_LSOLVE_FAIL, CV_RHSFUNC_FAIL
 *    RHSFUNC_RECVR,
 *    TRY_AGAIN
 *
 */

#define DO_ERROR_TEST +2
#define PREDICT_AGAIN +3

#define TRY_AGAIN      +5
#define FIRST_CALL     +6
#define PREV_CONV_FAIL +7
#define PREV_PROJ_FAIL +8
#define PREV_ERR_FAIL  +9

#define RHSFUNC_RECVR    +10
#define CONSTR_RECVR     +11
#define CONSTRFUNC_RECVR +12
#define PROJFUNC_RECVR   +13

#define QRHSFUNC_RECVR  +14
#define SRHSFUNC_RECVR  +15
#define QSRHSFUNC_RECVR +16

/* nonlinear solver constants
   NLS_MAXCOR  maximum no. of corrector iterations for the nonlinear solver
   CRDOWN      constant used in the estimation of the convergence rate (crate)
               of the iterates for the nonlinear equation
   RDIV        declare divergence if ratio del/delp > RDIV
*/
#define NLS_MAXCOR 3
#define CRDOWN     SUN_RCONST(0.3)
#define RDIV       SUN_RCONST(2.0)

/*
 * =================================================================
 *   F O R W A R D   P O I N T E R   R E F E R E N C E S
 * =================================================================
 */

typedef struct CVadjMemRec* CVadjMem;
typedef struct CVckpntMemRec* CVckpntMem;
typedef struct CVdtpntMemRec* CVdtpntMem;
typedef struct CVodeBMemRec* CVodeBMem;

/*
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types: struct CVodeMemRec, CVodeMem
 * -----------------------------------------------------------------
 * The type CVodeMem is type pointer to struct CVodeMemRec.
 * This structure contains fields to keep track of problem state.
 * -----------------------------------------------------------------
 */

typedef struct CVodeMemRec
{
  SUNContext cv_sunctx;

  void* python;

  sunrealtype cv_uround; /* machine unit roundoff */

  /*--------------------------
    Problem Specification Data
    --------------------------*/

  CVRhsFn cv_f;       /* y' = f(t,y(t))                                */
  void* cv_user_data; /* user pointer passed to f                      */
  int cv_lmm;         /* lmm = CV_ADAMS or CV_BDF                      */
  int cv_itol;        /* itol = CV_SS, CV_SV, CV_WF, CV_NN             */

  sunrealtype cv_reltol;  /* relative tolerance                            */
  sunrealtype cv_Sabstol; /* scalar absolute tolerance                     */
  N_Vector cv_Vabstol;    /* vector absolute tolerance                     */
  sunbooleantype cv_atolmin0; /* flag indicating that min(abstol) = 0          */
  sunbooleantype cv_user_efun; /* SUNTRUE if user sets efun                     */
  CVEwtFn cv_efun; /* function to set ewt                           */
  void* cv_e_data; /* user pointer passed to efun                   */

  sunbooleantype cv_constraintsSet; /* constraints vector present:
                                    do constraints calc                       */

  /*-----------------------
    Quadrature Related Data
    -----------------------*/

  sunbooleantype cv_quadr; /* SUNTRUE if integrating quadratures            */

  CVQuadRhsFn cv_fQ; /* q' = fQ(t, y(t))                              */

  sunbooleantype cv_errconQ; /* SUNTRUE if quadrs. are included in error test */

  int cv_itolQ;            /* itolQ = CV_SS or CV_SV                        */
  sunrealtype cv_reltolQ;  /* relative tolerance for quadratures            */
  sunrealtype cv_SabstolQ; /* scalar absolute tolerance for quadratures     */
  N_Vector cv_VabstolQ;    /* vector absolute tolerance for quadratures     */
  sunbooleantype cv_atolQmin0; /* flag indicating that min(abstolQ) = 0         */

  /*------------------------
    Sensitivity Related Data
    ------------------------*/

  sunbooleantype cv_sensi; /* SUNTRUE if computing sensitivities           */

  int cv_Ns; /* Number of sensitivities                      */

  int cv_ism; /* ism = SIMULTANEOUS or STAGGERED              */

  CVSensRhsFn cv_fS;      /* fS = (df/dy)*yS + (df/dp)                    */
  CVSensRhs1Fn cv_fS1;    /* fS1 = (df/dy)*yS_i + (df/dp)                 */
  void* cv_fS_data;       /* data pointer passed to fS                    */
  sunbooleantype cv_fSDQ; /* SUNTRUE if using internal DQ functions       */
  int cv_ifS;             /* ifS = ALLSENS or ONESENS                     */

  sunrealtype* cv_p;       /* parameters in f(t,y,p)                       */
  sunrealtype* cv_pbar;    /* scale factors for parameters                 */
  int* cv_plist;           /* list of sensitivities                        */
  int cv_DQtype;           /* central/forward finite differences           */
  sunrealtype cv_DQrhomax; /* cut-off value for separate/simultaneous FD   */

  sunbooleantype cv_errconS; /* SUNTRUE if yS are considered in err. control */

  int cv_itolS;
  sunrealtype cv_reltolS;   /* relative tolerance for sensitivities         */
  sunrealtype* cv_SabstolS; /* scalar absolute tolerances for sensi.        */
  N_Vector* cv_VabstolS;    /* vector absolute tolerances for sensi.        */
  sunbooleantype* cv_atolSmin0; /* flags indicating that min(abstolS[i]) = 0    */

  /*-----------------------------------
    Quadrature Sensitivity Related Data
    -----------------------------------*/

  sunbooleantype cv_quadr_sensi; /* SUNTRUE if computing sensitivities of quadrs. */

  CVQuadSensRhsFn cv_fQS;  /* fQS = (dfQ/dy)*yS + (dfQ/dp)                 */
  void* cv_fQS_data;       /* data pointer passed to fQS                   */
  sunbooleantype cv_fQSDQ; /* SUNTRUE if using internal DQ functions       */

  sunbooleantype cv_errconQS; /* SUNTRUE if yQS are considered in err. con.   */

  int cv_itolQS;
  sunrealtype cv_reltolQS;   /* relative tolerance for yQS                   */
  sunrealtype* cv_SabstolQS; /* scalar absolute tolerances for yQS           */
  N_Vector* cv_VabstolQS;    /* vector absolute tolerances for yQS           */
  sunbooleantype* cv_atolQSmin0; /* flags indicating that min(abstolQS[i]) = 0   */

  /*-----------------------
    Nordsieck History Array
    -----------------------*/

  N_Vector cv_zn[L_MAX]; /* Nordsieck array, of size N x (q+1).
                             zn[j] is a vector of length N (j=0,...,q)
                             zn[j] = [1/factorial(j)] * h^j *
                             (jth derivative of the interpolating polynomial) */

  /*-------------------
    Vectors of length N
    -------------------*/

  N_Vector cv_ewt;    /* error weight vector                                 */
  N_Vector cv_y;      /* y is used as temporary storage by the solver
                          The memory is provided by the user to CVode
                          where the vector is named yout.                     */
  N_Vector cv_acor;   /* In the context of the solution of the nonlinear
                          equation, acor = y_n(m) - y_n(0). On return, this
                          vector is scaled to give the estimated local error  */
  N_Vector cv_tempv;  /* temporary storage vector                            */
  N_Vector cv_ftemp;  /* temporary storage vector                            */
  N_Vector cv_vtemp1; /* temporary storage vector                            */
  N_Vector cv_vtemp2; /* temporary storage vector                            */
  N_Vector cv_vtemp3; /* temporary storage vector                            */

  N_Vector cv_constraints; /* vector of inequality constraint options         */

  /*--------------------------
    Quadrature Related Vectors
    --------------------------*/

  N_Vector cv_znQ[L_MAX]; /* Nordsieck arrays for quadratures             */
  N_Vector cv_ewtQ;       /* error weight vector for quadratures          */
  N_Vector cv_yQ;         /* Unlike y, yQ is not allocated by the user    */
  N_Vector cv_acorQ;      /* acorQ = yQ_n(m) - yQ_n(0)                    */
  N_Vector cv_tempvQ;     /* temporary storage vector (~ tempv)           */

  /*---------------------------
    Sensitivity Related Vectors
    ---------------------------*/

  N_Vector* cv_znS[L_MAX]; /* Nordsieck arrays for sensitivities           */
  N_Vector* cv_ewtS;       /* error weight vectors for sensitivities       */
  N_Vector* cv_yS;         /* yS=yS0 (allocated by the user)               */
  N_Vector* cv_acorS;      /* acorS = yS_n(m) - yS_n(0)                    */
  N_Vector* cv_tempvS;     /* temporary storage vector (~ tempv)           */
  N_Vector* cv_ftempS;     /* temporary storage vector (~ ftemp)           */

  sunbooleantype cv_stgr1alloc; /* Did we allocate ncfS1, ncfnS1, and nniS1?    */

  /*--------------------------------------
    Quadrature Sensitivity Related Vectors
    --------------------------------------*/

  N_Vector* cv_znQS[L_MAX]; /* Nordsieck arrays for quadr. sensitivities    */
  N_Vector* cv_ewtQS;       /* error weight vectors for sensitivities       */
  N_Vector* cv_yQS;         /* Unlike yS, yQS is not allocated by the user  */
  N_Vector* cv_acorQS;      /* acorQS = yQS_n(m) - yQS_n(0)                 */
  N_Vector* cv_tempvQS;     /* temporary storage vector (~ tempv)           */
  N_Vector cv_ftempQ;       /* temporary storage vector (~ ftemp)           */

  /*-----------------
    Tstop information
    -----------------*/

  sunbooleantype cv_tstopset;
  sunbooleantype cv_tstopinterp;
  sunrealtype cv_tstop;

  /*---------
    Step Data
    ---------*/

  int cv_q;      /* current order                               */
  int cv_qprime; /* order to be used on the next step
                                  qprime = q-1, q, or q+1                     */
  int cv_next_q; /* order to be used on the next step           */
  int cv_qwait;  /* number of internal steps to wait before
                                  considering a change in q                   */
  int cv_L;      /* L = q + 1                                   */

  sunrealtype cv_hin;      /* initial step size                           */
  sunrealtype cv_h;        /* current step size                           */
  sunrealtype cv_hprime;   /* step size to be used on the next step       */
  sunrealtype cv_next_h;   /* step size to be used on the next step       */
  sunrealtype cv_eta;      /* eta = hprime / h                            */
  sunrealtype cv_hscale;   /* value of h used in zn                       */
  sunrealtype cv_tn;       /* current internal value of t                 */
  sunrealtype cv_tretlast; /* last value of t returned by CVode           */

  sunrealtype cv_tau[L_MAX + 1];    /* array of previous q+1 successful step
                                  sizes indexed from 1 to q+1                 */
  sunrealtype cv_tq[NUM_TESTS + 1]; /* array of test quantities indexed from
                                  1 to NUM_TESTS(=5)                          */
  sunrealtype cv_l[L_MAX]; /* coefficients of l(x) (degree q poly)        */

  sunrealtype cv_rl1;    /* the scalar 1/l[1]                           */
  sunrealtype cv_gamma;  /* gamma = h * rl1                             */
  sunrealtype cv_gammap; /* gamma at the last setup call                */
  sunrealtype cv_gamrat; /* gamma / gammap                              */

  sunrealtype cv_crate;        /* estimated corrector convergence rate        */
  sunrealtype cv_crateS;       /* estimated corrector convergence rate (Stgr) */
  sunrealtype cv_delp;         /* norm of previous nonlinear solver update    */
  sunrealtype cv_acnrm;        /* | acor |                                    */
  sunbooleantype cv_acnrmcur;  /* is | acor | current?                        */
  sunrealtype cv_acnrmQ;       /* | acorQ |                                   */
  sunrealtype cv_acnrmS;       /* | acorS |                                   */
  sunbooleantype cv_acnrmScur; /* is | acorS | current?                       */
  sunrealtype cv_acnrmQS;      /* | acorQS |                                  */
  sunrealtype cv_nlscoef;      /* coefficient in nonlinear convergence test   */
  int* cv_ncfS1;               /* Array of Ns local counters for conv.
                                * failures (used in CVStep for STAGGERED1)    */

  /*------
    Limits
    ------*/

  int cv_qmax;        /* q <= qmax                                          */
  long int cv_mxstep; /* maximum number of internal steps for one user call */
  int cv_mxhnil;      /* maximum number of warning messages issued to the
                             user that t + h == t for the next internal step    */
  int cv_maxnef;      /* maximum number of error test failures              */
  int cv_maxncf;      /* maximum number of nonlinear convergence failures   */

  sunrealtype cv_hmin; /* |h| >= hmin                                        */
  sunrealtype cv_hmax_inv; /* |h| <= 1/hmax_inv                                  */
  sunrealtype cv_etamax; /* eta <= etamax                                      */
  sunrealtype cv_eta_min_fx; /* eta_min_fx < eta < eta_max_fx keep the current h   */
  sunrealtype cv_eta_max_fx;
  sunrealtype cv_eta_max_fs; /* eta <= eta_max_fs on the first step                */
  sunrealtype cv_eta_max_es; /* eta <= eta_max_es on early steps                   */
  sunrealtype cv_eta_max_gs; /* eta <= eta_max_gs on a general step                */
  sunrealtype cv_eta_min; /* eta >= eta_min on a general step                   */
  sunrealtype cv_eta_min_ef; /* eta >= eta_min_ef after an error test failure      */
  sunrealtype cv_eta_max_ef; /* eta on multiple (>= small_nef) error test failures */
  sunrealtype cv_eta_cf; /* eta on a nonlinear solver convergence failure      */

  long int cv_small_nst; /* nst <= small_nst use eta_max_es */
  int cv_small_nef;      /* nef >= small_nef use eta_max_ef */

  /*--------
    Counters
    --------*/

  long int cv_nst; /* number of internal steps taken                  */

  long int cv_nfe;   /* number of f calls                               */
  long int cv_nfQe;  /* number of fQ calls                              */
  long int cv_nfSe;  /* number of fS calls                              */
  long int cv_nfeS;  /* number of f calls from sensi DQ                 */
  long int cv_nfQSe; /* number of fQS calls                             */
  long int cv_nfQeS; /* number of fQ calls from sensi DQ                */

  long int cv_ncfn;    /* number of corrector convergence failures        */
  long int cv_ncfnS;   /* number of total sensi. corr. conv. failures     */
  long int* cv_ncfnS1; /* number of sensi. corrector conv. failures       */

  long int cv_nni;    /* number of nonlinear iterations performed        */
  long int cv_nniS;   /* number of total sensi. nonlinear iterations     */
  long int* cv_nniS1; /* number of sensi. nonlinear iterations           */

  long int cv_nnf;    /* number of nonlinear convergence fails           */
  long int cv_nnfS;   /* number of total sensi. nonlinear conv. fails    */
  long int* cv_nnfS1; /* number of sensi. nonlinear conv. fails          */

  long int cv_netf;   /* number of error test failures                   */
  long int cv_netfQ;  /* number of quadr. error test failures            */
  long int cv_netfS;  /* number of sensi. error test failures            */
  long int cv_netfQS; /* number of quadr. sensi. error test failures     */

  long int cv_nsetups;  /* number of setup calls                           */
  long int cv_nsetupsS; /* number of setup calls due to sensitivities      */

  int cv_nhnil; /* number of messages issued to the user that
                              t + h == t for the next iternal step            */

  /*----------------
    Step size ratios
    ----------------*/

  sunrealtype cv_etaqm1; /* ratio of new to old h for order q-1             */
  sunrealtype cv_etaq;   /* ratio of new to old h for order q               */
  sunrealtype cv_etaqp1; /* ratio of new to old h for order q+1             */

  /*------------------
    Space requirements
    ------------------*/

  sunindextype cv_lrw1; /* no. of sunrealtype words in 1 N_Vector y           */
  sunindextype cv_liw1; /* no. of integer words in 1 N_Vector y            */
  sunindextype cv_lrw1Q; /* no. of sunrealtype words in 1 N_Vector yQ          */
  sunindextype cv_liw1Q; /* no. of integer words in 1 N_Vector yQ           */
  long int cv_lrw; /* no. of sunrealtype words in CVODE work vectors     */
  long int cv_liw; /* no. of integer words in CVODE work vectors      */

  /*---------------------
    Nonlinear Solver Data
    ---------------------*/

  SUNNonlinearSolver NLS; /* nonlinear solver object                   */
  sunbooleantype ownNLS;  /* flag indicating NLS ownership             */

  SUNNonlinearSolver NLSsim; /* NLS object for the simultaneous corrector */
  sunbooleantype ownNLSsim;  /* flag indicating NLS ownership             */

  SUNNonlinearSolver NLSstg; /* NLS object for the staggered corrector */
  sunbooleantype ownNLSstg;  /* flag indicating NLS ownership          */

  SUNNonlinearSolver NLSstg1; /* NLS object for the staggered1 corrector */
  sunbooleantype ownNLSstg1;  /* flag indicating NLS ownership           */
  int sens_solve_idx;         /* index of the current staggered1 solve   */
  long int nnip;              /* previous total number of iterations     */

  sunbooleantype sens_solve; /* flag indicating if the current solve is a
                                  staggered or staggered1 sensitivity solve */
  CVRhsFn nls_f;             /* f(t,y(t)) used in the nonlinear solver    */
  int convfail;              /* flag to indicate when a Jacobian update may
                                  be needed */

  /* The following vectors are NVector wrappers for use with the simultaneous
     and staggered corrector methods:

       Simultaneous: zn0Sim  = [cv_zn[0], cv_znS[0]]
                     ycorSim = [cv_acor,  cv_acorS]
                     ewtSim  = [cv_ewt,   cv_ewtS]

       Staggered: zn0Stg  = cv_znS[0]
                  ycorStg = cv_acorS
                  ewtStg  = cv_ewtS
  */
  N_Vector zn0Sim, ycorSim, ewtSim;
  N_Vector zn0Stg, ycorStg, ewtStg;

  /* flags indicating if vector wrappers for the simultaneous and staggered
     correctors have been allocated */
  sunbooleantype simMallocDone;
  sunbooleantype stgMallocDone;

  /*------------------
    Linear Solver Data
    ------------------*/

  /* Linear Solver functions to be called */

  int (*cv_linit)(struct CVodeMemRec* cv_mem);

  int (*cv_lreinit)(struct CVodeMemRec* cv_mem);

  int (*cv_lsetup)(struct CVodeMemRec* cv_mem, int convfail, N_Vector ypred,
                   N_Vector fpred, sunbooleantype* jcurPtr, N_Vector vtemp1,
                   N_Vector vtemp2, N_Vector vtemp3);

  int (*cv_lsolve)(struct CVodeMemRec* cv_mem, N_Vector b, N_Vector weight,
                   N_Vector ycur, N_Vector fcur);

  int (*cv_lfree)(struct CVodeMemRec* cv_mem);

  /* Linear Solver specific memory */

  void* cv_lmem;                /* linear solver interface memory structure */
  long int cv_msbp;             /* max number of steps between lsetip calls */
  sunrealtype cv_dgmax_lsetup;  /* gamma ratio threshold to signal for a linear
                              * solver setup */
  sunbooleantype cv_forceSetup; /* flag to request a call to the setup routine */

  /*------------
    Saved Values
    ------------*/

  int cv_qu;                /* last successful q value used                */
  long int cv_nstlp;        /* step number of last setup call              */
  sunrealtype cv_h0u;       /* actual initial stepsize                     */
  sunrealtype cv_hu;        /* last successful h value used                */
  sunrealtype cv_saved_tq5; /* saved value of tq[5]                        */
  sunbooleantype cv_jcur;   /* is Jacobian info for linear solver current? */
  int cv_convfail;          /* flag storing previous solver failure mode   */
  sunrealtype cv_tolsf;     /* tolerance scale factor                      */
  int cv_qmax_alloc;        /* value of qmax used when allocating mem      */
  int cv_qmax_allocQ;       /* qmax used when allocating quad. mem         */
  int cv_qmax_allocS;       /* qmax used when allocating sensi. mem        */
  int cv_qmax_allocQS;      /* qmax used when allocating quad. sensi. mem  */
  int cv_indx_acor;         /* index of the zn vector with saved acor      */

  /*--------------------------------------------------------------------
    Flags turned ON by CVodeInit, CVodeSensMalloc, and CVodeQuadMalloc
    and read by CVodeReInit, CVodeSensReInit, and CVodeQuadReInit
    --------------------------------------------------------------------*/

  sunbooleantype cv_VabstolMallocDone;
  sunbooleantype cv_MallocDone;
  sunbooleantype cv_constraintsMallocDone;

  sunbooleantype cv_VabstolQMallocDone;
  sunbooleantype cv_QuadMallocDone;

  sunbooleantype cv_VabstolSMallocDone;
  sunbooleantype cv_SabstolSMallocDone;
  sunbooleantype cv_SensMallocDone;

  sunbooleantype cv_VabstolQSMallocDone;
  sunbooleantype cv_SabstolQSMallocDone;
  sunbooleantype cv_QuadSensMallocDone;

  /*-------------------------------------------
    User access function
    -------------------------------------------*/
  CVMonitorFn cv_monitorfun;    /* func called with CVODE mem and user data  */
  long int cv_monitor_interval; /* step interval to call cv_monitorfun       */

  /*-------------------------
    Stability Limit Detection
    -------------------------*/

  sunbooleantype cv_sldeton;  /* is Stability Limit Detection on?             */
  sunrealtype cv_ssdat[6][4]; /* scaled data array for STALD                  */
  int cv_nscon;               /* counter for STALD method                     */
  long int cv_nor;            /* counter for number of order reductions       */

  /*----------------
    Rootfinding Data
    ----------------*/

  CVRootFn cv_gfun;      /* function g for roots sought                     */
  int cv_nrtfn;          /* number of components of g                       */
  int* cv_iroots;        /* array for root information                      */
  int* cv_rootdir;       /* array specifying direction of zero-crossing     */
  sunrealtype cv_tlo;    /* nearest endpoint of interval in root search     */
  sunrealtype cv_thi;    /* farthest endpoint of interval in root search    */
  sunrealtype cv_trout;  /* t value returned by rootfinding routine         */
  sunrealtype* cv_glo;   /* saved array of g values at t = tlo              */
  sunrealtype* cv_ghi;   /* saved array of g values at t = thi              */
  sunrealtype* cv_grout; /* array of g values at t = trout                  */
  sunrealtype cv_ttol;   /* tolerance on root location trout                */
  int cv_irfnd;          /* flag showing whether last step had a root       */
  long int cv_nge;       /* counter for g evaluations                       */
  sunbooleantype* cv_gactive; /* array with active/inactive event functions      */
  int cv_mxgnull; /* number of warning messages about possible g==0  */

  /*---------------
    Projection Data
    ---------------*/

  CVodeProjMem proj_mem;       /* projection memory structure               */
  sunbooleantype proj_enabled; /* flag indicating if projection is enabled  */
  sunbooleantype proj_applied; /* flag indicating if projection was applied */
  sunrealtype proj_p[L_MAX];   /* coefficients of p(x) (degree q poly)      */

  /*-----------------------
    Fused Vector Operations
    -----------------------*/

  sunrealtype* cv_cvals; /* array of scalars */
  N_Vector* cv_Xvecs;    /* array of vectors */
  N_Vector* cv_Zvecs;    /* array of vectors */

  /*----------------
    Resizing History
    ----------------*/

  sunbooleantype first_step_after_resize; /* Flag to signal a resize happened */

  /*------------------------
    Adjoint sensitivity data
    ------------------------*/

  sunbooleantype cv_adj; /* SUNTRUE if performing ASA                */

  struct CVadjMemRec* cv_adj_mem; /* Pointer to adjoint memory structure      */

  sunbooleantype cv_adjMallocDone;

}* CVodeMem;

/*
 * =================================================================
 *   A D J O I N T   M O D U L E    M E M O R Y    B L O C K
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : struct CVckpntMemRec, CVckpntMem
 * -----------------------------------------------------------------
 * The type CVckpntMem is type pointer to struct CVckpntMemRec.
 * This structure contains fields to store all information at a
 * check point that is needed to 'hot' start cvodes.
 * -----------------------------------------------------------------
 */

struct CVckpntMemRec
{
  /* Integration limits */
  sunrealtype ck_t0;
  sunrealtype ck_t1;

  /* Nordsieck History Array */
  N_Vector ck_zn[L_MAX];

  /* Do we need to carry quadratures? */
  sunbooleantype ck_quadr;

  /* Nordsieck History Array for quadratures */
  N_Vector ck_znQ[L_MAX];

  /* Do we need to carry sensitivities? */
  sunbooleantype ck_sensi;

  /* number of sensitivities */
  int ck_Ns;

  /* Nordsieck History Array for sensitivities */
  N_Vector* ck_znS[L_MAX];

  /* Do we need to carry quadrature sensitivities? */
  sunbooleantype ck_quadr_sensi;

  /* Nordsieck History Array for quadrature sensitivities */
  N_Vector* ck_znQS[L_MAX];

  /* Was ck_zn[qmax] allocated?
     ck_zqm = 0    - no
     ck_zqm = qmax - yes      */
  int ck_zqm;

  /* Step data */
  long int ck_nst;
  sunrealtype ck_tretlast;
  int ck_q;
  int ck_qprime;
  int ck_qwait;
  int ck_L;
  sunrealtype ck_gammap;
  sunrealtype ck_h;
  sunrealtype ck_hprime;
  sunrealtype ck_hscale;
  sunrealtype ck_eta;
  sunrealtype ck_etamax;
  sunrealtype ck_tau[L_MAX + 1];
  sunrealtype ck_tq[NUM_TESTS + 1];
  sunrealtype ck_l[L_MAX];

  /* Saved values */
  sunrealtype ck_saved_tq5;

  /* Pointer to next structure in list */
  struct CVckpntMemRec* ck_next;
};

/*
 * -----------------------------------------------------------------
 * Types for functions provided by an interpolation module
 * -----------------------------------------------------------------
 * cvaIMMallocFn: Type for a function that initializes the content
 *                field of the structures in the dt array
 * cvaIMFreeFn:   Type for a function that deallocates the content
 *                field of the structures in the dt array
 * cvaIMGetYFn:   Type for a function that returns the
 *                interpolated forward solution.
 * cvaIMStorePnt: Type for a function that stores a new
 *                point in the structure d
 * -----------------------------------------------------------------
 */

typedef sunbooleantype (*cvaIMMallocFn)(CVodeMem cv_mem);
typedef void (*cvaIMFreeFn)(CVodeMem cv_mem);
typedef int (*cvaIMGetYFn)(CVodeMem cv_mem, sunrealtype t, N_Vector y,
                           N_Vector* yS);
typedef int (*cvaIMStorePntFn)(CVodeMem cv_mem, CVdtpntMem d);

/*
 * -----------------------------------------------------------------
 * Type : struct CVdtpntMemRec
 * -----------------------------------------------------------------
 * This structure contains fields to store all information at a
 * data point that is needed to interpolate solution of forward
 * simulations. Its content field depends on IMtype.
 * -----------------------------------------------------------------
 */

struct CVdtpntMemRec
{
  sunrealtype t; /* time */
  void* content; /* IMtype-dependent content */
};

/* Data for cubic Hermite interpolation */
typedef struct CVhermiteDataMemRec
{
  N_Vector y;
  N_Vector yd;
  N_Vector* yS;
  N_Vector* ySd;
}* CVhermiteDataMem;

/* Data for polynomial interpolation */
typedef struct CVpolynomialDataMemRec
{
  N_Vector y;
  N_Vector* yS;
  int order;
}* CVpolynomialDataMem;

/*
 * -----------------------------------------------------------------
 * Type : struct CVodeBMemRec
 * -----------------------------------------------------------------
 * The type CVodeBMem is a pointer to a structure which stores all
 * information for ONE backward problem.
 * The CVadjMem structure contains a linked list of CVodeBMem pointers
 * -----------------------------------------------------------------
 */

struct CVodeBMemRec
{
  /* Index of this backward problem */
  int cv_index;

  /* Time at which the backward problem is initialized */
  sunrealtype cv_t0;

  /* CVODES memory for this backward problem */
  CVodeMem cv_mem;

  /* Flags to indicate that this backward problem's RHS or quad RHS
   * require forward sensitivities */
  sunbooleantype cv_f_withSensi;
  sunbooleantype cv_fQ_withSensi;

  /* Right hand side function for backward run */
  CVRhsFnB cv_f;
  CVRhsFnBS cv_fs;

  /* Right hand side quadrature function for backward run */
  CVQuadRhsFnB cv_fQ;
  CVQuadRhsFnBS cv_fQs;

  /* User user_data */
  void* cv_user_data;

  /* Memory block for a linear solver's interface to CVODEA */
  void* cv_lmem;

  /* Function to free any memory allocated by the linear solver */
  int (*cv_lfree)(CVodeBMem cvB_mem);

  /* Memory block for a preconditioner's module interface to CVODEA */
  void* cv_pmem;

  /* Function to free any memory allocated by the preconditioner module */
  int (*cv_pfree)(CVodeBMem cvB_mem);

  /* Time at which to extract solution / quadratures */
  sunrealtype cv_tout;

  /* Workspace Nvector */
  N_Vector cv_y;

  /* Pointer to next structure in list */
  struct CVodeBMemRec* cv_next;
};

/*
 * -----------------------------------------------------------------
 * Type : struct CVadjMemRec
 * -----------------------------------------------------------------
 * The type CVadjMem is type pointer to struct CVadjMemRec.
 * This structure contains fields to store all information
 * necessary for adjoint sensitivity analysis.
 * -----------------------------------------------------------------
 */

struct CVadjMemRec
{
  /* --------------------
   * Forward problem data
   * -------------------- */

  /* Integration interval */
  sunrealtype ca_tinitial, ca_tfinal;

  /* Flag for first call to CVodeF */
  sunbooleantype ca_firstCVodeFcall;

  /* Flag if CVodeF was called with TSTOP */
  sunbooleantype ca_tstopCVodeFcall;
  sunrealtype ca_tstopCVodeF;

  /* Flag if CVodeF was called in CV_NORMAL_MODE and encountered a
     root after tout */
  sunbooleantype ca_rootret;
  sunrealtype ca_troot;

  /* ----------------------
   * Backward problems data
   * ---------------------- */

  /* Storage for backward problems */
  struct CVodeBMemRec* cvB_mem;

  /* Number of backward problems */
  int ca_nbckpbs;

  /* Address of current backward problem */
  struct CVodeBMemRec* ca_bckpbCrt;

  /* Flag for first call to CVodeB */
  sunbooleantype ca_firstCVodeBcall;

  /* ----------------
   * Check point data
   * ---------------- */

  /* Storage for check point information */
  struct CVckpntMemRec* ck_mem;

  /* Number of check points */
  int ca_nckpnts;

  /* address of the check point structure for which data is available */
  struct CVckpntMemRec* ca_ckpntData;

  /* ------------------
   * Interpolation data
   * ------------------ */

  /* Number of steps between 2 check points */
  long int ca_nsteps;

  /* Last index used in CVAfindIndex */
  long int ca_ilast;

  /* Storage for data from forward runs */
  struct CVdtpntMemRec** dt_mem;

  /* Actual number of data points in dt_mem (typically np=nsteps+1) */
  long int ca_np;

  /* Interpolation type */
  int ca_IMtype;

  /* Functions set by the interpolation module */
  cvaIMMallocFn ca_IMmalloc;
  cvaIMFreeFn ca_IMfree;
  cvaIMStorePntFn ca_IMstore; /* store a new interpolation point */
  cvaIMGetYFn ca_IMget;       /* interpolate forward solution    */

  /* Flags controlling the interpolation module */
  sunbooleantype ca_IMmallocDone;  /* IM initialized? */
  sunbooleantype ca_IMnewData;     /* new data available in dt_mem?*/
  sunbooleantype ca_IMstoreSensi;  /* store sensitivities? */
  sunbooleantype ca_IMinterpSensi; /* interpolate sensitivities? */

  /* Workspace for the interpolation module */
  N_Vector ca_Y[L_MAX];   /* pointers to zn[i] */
  N_Vector* ca_YS[L_MAX]; /* pointers to znS[i] */
  sunrealtype ca_T[L_MAX];

  /* -------------------------------
   * Workspace for wrapper functions
   * ------------------------------- */

  N_Vector ca_ytmp;

  N_Vector* ca_yStmp;
};

/*
 * =================================================================
 *     I N T E R F A C E   T O    L I N E A R   S O L V E R S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Communication between CVODE and a CVODE Linear Solver
 * -----------------------------------------------------------------
 * convfail (input to cv_lsetup)
 *
 * CV_NO_FAILURES : Either this is the first cv_setup call for this
 *                  step, or the local error test failed on the
 *                  previous attempt at this step (but the nonlinear
 *                  solver iteration converged).
 *
 * CV_FAIL_BAD_J  : This value is passed to cv_lsetup if
 *
 *                  (a) The previous nonlinear solver corrector iteration
 *                      did not converge and the linear solver's
 *                      setup routine indicated that its Jacobian-
 *                      related data is not current
 *                                   or
 *                  (b) During the previous nonlinear solver corrector
 *                      iteration, the linear solver's solve routine
 *                      failed in a recoverable manner and the
 *                      linear solver's setup routine indicated that
 *                      its Jacobian-related data is not current.
 *
 * CV_FAIL_OTHER  : During the current internal step try, the
 *                  previous nonlinear solver iteration failed to converge
 *                  even though the linear solver was using current
 *                  Jacobian-related data.
 * -----------------------------------------------------------------
 */

/* Constants for convfail (input to cv_lsetup) */

#define CV_NO_FAILURES 0
#define CV_FAIL_BAD_J  1
#define CV_FAIL_OTHER  2

/*
 * -----------------------------------------------------------------
 * int (*cv_linit)(CVodeMem cv_mem);
 * -----------------------------------------------------------------
 * The purpose of cv_linit is to complete initializations for a
 * specific linear solver, such as counters and statistics.
 * An LInitFn should return 0 if it has successfully initialized the
 * CVODE linear solver and a negative value otherwise.
 * If an error does occur, an appropriate message should be sent to
 * the error handler function.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cv_lsetup)(CVodeMem cv_mem, int convfail, N_Vector ypred,
 *                 N_Vector fpred, sunbooleantype *jcurPtr,
 *                 N_Vector vtemp1, N_Vector vtemp2,
 *                 N_Vector vtemp3);
 * -----------------------------------------------------------------
 * The job of cv_lsetup is to prepare the linear solver for
 * subsequent calls to cv_lsolve. It may recompute Jacobian-
 * related data is it deems necessary. Its parameters are as
 * follows:
 *
 * cv_mem - problem memory pointer of type CVodeMem. See the
 *          typedef earlier in this file.
 *
 * convfail - a flag to indicate any problem that occurred during
 *            the solution of the nonlinear equation on the
 *            current time step for which the linear solver is
 *            being used. This flag can be used to help decide
 *            whether the Jacobian data kept by a CVODE linear
 *            solver needs to be updated or not.
 *            Its possible values have been documented above.
 *
 * ypred - the predicted y vector for the current CVODE internal
 *         step.
 *
 * fpred - f(tn, ypred).
 *
 * jcurPtr - a pointer to a boolean to be filled in by cv_lsetup.
 *           The function should set *jcurPtr=SUNTRUE if its Jacobian
 *           data is current after the call and should set
 *           *jcurPtr=SUNFALSE if its Jacobian data is not current.
 *           Note: If cv_lsetup calls for re-evaluation of
 *           Jacobian data (based on convfail and CVODE state
 *           data), it should return *jcurPtr=SUNTRUE always;
 *           otherwise an infinite loop can result.
 *
 * vtemp1 - temporary N_Vector provided for use by cv_lsetup.
 *
 * vtemp3 - temporary N_Vector provided for use by cv_lsetup.
 *
 * vtemp3 - temporary N_Vector provided for use by cv_lsetup.
 *
 * The cv_lsetup routine should return 0 if successful, a positive
 * value for a recoverable error, and a negative value for an
 * unrecoverable error.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cv_lsolve)(CVodeMem cv_mem, N_Vector b, N_Vector weight,
 *                  N_Vector ycur, N_Vector fcur);
 * -----------------------------------------------------------------
 * cv_lsolve must solve the linear equation P x = b, where
 * P is some approximation to (I - gamma J), J = (df/dy)(tn,ycur)
 * and the RHS vector b is input. The N-vector ycur contains
 * the solver's current approximation to y(tn) and the vector
 * fcur contains the N_Vector f(tn,ycur). The solution is to be
 * returned in the vector b. cv_lsolve returns a positive value
 * for a recoverable error and a negative value for an
 * unrecoverable error. Success is indicated by a 0 return value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cv_lfree)(CVodeMem cv_mem);
 * -----------------------------------------------------------------
 * cv_lfree should free up any memory allocated by the linear
 * solver. This routine is called once a problem has been
 * completed and the linear solver is no longer needed.  It should
 * return 0 upon success, nonzero on failure.
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 *    I N T E R N A L   F U N C T I O N S
 * =================================================================
 */

/* Norm functions */

sunrealtype cvSensNorm(CVodeMem cv_mem, N_Vector* xS, N_Vector* wS);

sunrealtype cvSensUpdateNorm(CVodeMem cv_mem, sunrealtype old_nrm, N_Vector* xS,
                             N_Vector* wS);

/* Prototype of internal ewtSet function */

int cvEwtSet(N_Vector ycur, N_Vector weight, void* data);

/* High level error handler */

void cvProcessError(CVodeMem cv_mem, int error_code, int line, const char* func,
                    const char* file, const char* msgfmt, ...);

/* Prototype of internal ErrHandler function */

void cvErrHandler(int error_code, const char* module, const char* function,
                  char* msg, void* data);

/* Nonlinear solver initialization */

int cvNlsInit(CVodeMem cv_mem);
int cvNlsInitSensSim(CVodeMem cv_mem);
int cvNlsInitSensStg(CVodeMem cv_mem);
int cvNlsInitSensStg1(CVodeMem cv_mem);

/* Projection functions */

int cvDoProjection(CVodeMem cv_mem, int* nflagPtr, sunrealtype saved_t,
                   int* npfPtr);
int cvProjInit(CVodeProjMem proj_mem);
int cvProjFree(CVodeProjMem* proj_mem);

/* Restore tn and undo prediction to reattempt a step */

void cvRestore(CVodeMem cv_mem, sunrealtype saved_t);

/* Reset h and rescale history array to prepare for a step */

void cvRescale(CVodeMem cv_mem);

/* Prototypes for internal sensitivity rhs wrappers */

int cvSensRhsWrapper(CVodeMem cv_mem, sunrealtype time, N_Vector ycur,
                     N_Vector fcur, N_Vector* yScur, N_Vector* fScur,
                     N_Vector temp1, N_Vector temp2);

int cvSensRhs1Wrapper(CVodeMem cv_mem, sunrealtype time, N_Vector ycur,
                      N_Vector fcur, int is, N_Vector yScur, N_Vector fScur,
                      N_Vector temp1, N_Vector temp2);

/* Prototypes for internal sensitivity rhs DQ functions */

int cvSensRhsInternalDQ(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                        N_Vector* yS, N_Vector* ySdot, void* fS_data,
                        N_Vector tempv, N_Vector ftemp);

int cvSensRhs1InternalDQ(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                         int is, N_Vector yS, N_Vector ySdot, void* fS_data,
                         N_Vector tempv, N_Vector ftemp);

/* Function to destroy function table allocated by the Python binding code */

#if defined(SUNDIALS_ENABLE_PYTHON)
void cvode_user_supplied_fn_table_destroy(void* ptr);
#endif

/*
 * =================================================================
 *    E R R O R    M E S S A G E S
 * =================================================================
 */

#define MSG_TIME   "t = " SUN_FORMAT_G
#define MSG_TIME_H "t = " SUN_FORMAT_G " and h = " SUN_FORMAT_G
#define MSG_TIME_INT                                                \
  "t = " SUN_FORMAT_G " is not between tcur - hold = " SUN_FORMAT_G \
  " and tcur = " SUN_FORMAT_G
#define MSG_TIME_TOUT  "tout = " SUN_FORMAT_G
#define MSG_TIME_TSTOP "tstop = " SUN_FORMAT_G

/* Initialization and I/O error messages */

#define MSGCV_NO_MEM     "cvode_mem = NULL illegal."
#define MSGCV_CVMEM_FAIL "Allocation of cvode_mem failed."
#define MSGCV_MEM_FAIL   "A memory request failed."
#define MSGCV_BAD_LMM \
  "Illegal value for lmm. The legal values are CV_ADAMS and CV_BDF."
#define MSGCV_NULL_SUNCTX "sunctx = NULL illegal."
#define MSGCV_NO_MALLOC   "Attempt to call before CVodeInit."
#define MSGCV_NEG_MAXORD  "maxord <= 0 illegal."
#define MSGCV_BAD_MAXORD  "Illegal attempt to increase maximum method order."
#define MSGCV_SET_SLDET \
  "Attempt to use stability limit detection with the CV_ADAMS method illegal."
#define MSGCV_NEG_HMIN       "hmin < 0 illegal."
#define MSGCV_NEG_HMAX       "hmax < 0 illegal."
#define MSGCV_BAD_HMIN_HMAX  "Inconsistent step size limits: hmin > hmax."
#define MSGCV_BAD_RELTOL     "reltol < 0 illegal."
#define MSGCV_BAD_ABSTOL     "abstol has negative component(s) (illegal)."
#define MSGCV_NULL_ABSTOL    "abstol = NULL illegal."
#define MSGCV_NULL_Y0        "y0 = NULL illegal."
#define MSGCV_Y0_FAIL_CONSTR "y0 fails to satisfy constraints."
#define MSGCV_BAD_ISM_CONSTR                                                \
  "Constraints can not be enforced while forward sensitivity is used with " \
  "simultaneous method"
#define MSGCV_NULL_F        "f = NULL illegal."
#define MSGCV_NULL_G        "g = NULL illegal."
#define MSGCV_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGCV_BAD_CONSTR    "Illegal values in constraints vector."
#define MSGCV_BAD_K         "Illegal value for k."
#define MSGCV_NULL_DKY      "dky = NULL illegal."
#define MSGCV_BAD_T         "Illegal value for t." MSG_TIME_INT
#define MSGCV_NO_ROOT       "Rootfinding was not initialized."
#define MSGCV_NLS_INIT_FAIL "The nonlinear solver's init routine failed."

#define MSGCV_NO_QUAD "Quadrature integration not activated."
#define MSGCV_BAD_ITOLQ \
  "Illegal value for itolQ. The legal values are CV_SS and CV_SV."
#define MSGCV_NULL_ABSTOLQ "abstolQ = NULL illegal."
#define MSGCV_BAD_RELTOLQ  "reltolQ < 0 illegal."
#define MSGCV_BAD_ABSTOLQ  "abstolQ has negative component(s) (illegal)."

#define MSGCV_SENSINIT_2 "Sensitivity analysis already initialized."
#define MSGCV_NO_SENSI   "Forward sensitivity analysis not activated."
#define MSGCV_BAD_ITOLS \
  "Illegal value for itolS. The legal values are CV_SS, CV_SV, and CV_EE."
#define MSGCV_NULL_ABSTOLS "abstolS = NULL illegal."
#define MSGCV_BAD_RELTOLS  "reltolS < 0 illegal."
#define MSGCV_BAD_ABSTOLS  "abstolS has negative component(s) (illegal)."
#define MSGCV_BAD_PBAR     "pbar has zero component(s) (illegal)."
#define MSGCV_BAD_PLIST    "plist has negative component(s) (illegal)."
#define MSGCV_BAD_NS       "NS <= 0 illegal."
#define MSGCV_NULL_YS0     "yS0 = NULL illegal."
#define MSGCV_BAD_ISM                                                       \
  "Illegal value for ism. Legal values are: CV_SIMULTANEOUS, CV_STAGGERED " \
  "and CV_STAGGERED1."
#define MSGCV_BAD_IFS \
  "Illegal value for ifS. Legal values are: CV_ALLSENS and CV_ONESENS."
#define MSGCV_BAD_ISM_IFS "Illegal ism = CV_STAGGERED1 for CVodeSensInit."
#define MSGCV_BAD_IS      "Illegal value for is."
#define MSGCV_NULL_DKYA   "dkyA = NULL illegal."
#define MSGCV_BAD_DQTYPE \
  "Illegal value for DQtype. Legal values are: CV_CENTERED and CV_FORWARD."
#define MSGCV_BAD_DQRHO "DQrhomax < 0 illegal."

#define MSGCV_BAD_ITOLQS \
  "Illegal value for itolQS. The legal values are CV_SS, CV_SV, and CV_EE."
#define MSGCV_NULL_ABSTOLQS "abstolQS = NULL illegal."
#define MSGCV_BAD_RELTOLQS  "reltolQS < 0 illegal."
#define MSGCV_BAD_ABSTOLQS  "abstolQS has negative component(s) (illegal)."
#define MSGCV_NO_QUADSENSI \
  "Forward sensitivity analysis for quadrature variables not activated."
#define MSGCV_NULL_YQS0 "yQS0 = NULL illegal."

/* CVode Error Messages */

#define MSGCV_NO_TOL      "No integration tolerances have been specified."
#define MSGCV_LSOLVE_NULL "The linear solver's solve routine is NULL."
#define MSGCV_YOUT_NULL   "yout = NULL illegal."
#define MSGCV_TRET_NULL   "tret = NULL illegal."
#define MSGCV_BAD_EWT     "Initial ewt has component(s) equal to zero (illegal)."
#define MSGCV_EWT_NOW_BAD "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSGCV_BAD_ITASK   "Illegal value for itask."
#define MSGCV_BAD_H0      "h0 and tout - t0 inconsistent."
#define MSGCV_BAD_TOUT                      \
  "Trouble interpolating at " MSG_TIME_TOUT \
  ". tout too far back in direction of integration"
#define MSGCV_EWT_FAIL "The user-provide EwtSet function failed."
#define MSGCV_EWT_NOW_FAIL \
  "At " MSG_TIME ", the user-provide EwtSet function failed."
#define MSGCV_LINIT_FAIL "The linear solver's init routine failed."
#define MSGCV_HNIL_DONE                                                    \
  "The above warning has been issued mxhnil times and will not be issued " \
  "again for this problem."
#define MSGCV_TOO_CLOSE "tout too close to t0 to start integration."
#define MSGCV_MAX_STEPS \
  "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSGCV_TOO_MUCH_ACC "At " MSG_TIME ", too much accuracy requested."
#define MSGCV_HNIL                                                         \
  "Internal " MSG_TIME_H " are such that t + h = t on the next step. The " \
  "solver will continue anyway."
#define MSGCV_ERR_FAILS \
  "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSGCV_CONV_FAILS \
  "At " MSG_TIME_H       \
  ", the corrector convergence test failed repeatedly or with |h| = hmin."
#define MSGCV_SETUP_FAILED \
  "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSGCV_SOLVE_FAILED \
  "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSGCV_FAILED_CONSTR \
  "At " MSG_TIME ", unable to satisfy inequality constraints."
#define MSGCV_RHSFUNC_FAILED \
  "At " MSG_TIME             \
  ", the right-hand side routine failed in an unrecoverable manner."
#define MSGCV_RHSFUNC_UNREC                                                   \
  "At " MSG_TIME ", the right-hand side failed in a recoverable manner, but " \
  "no recovery is possible."
#define MSGCV_RHSFUNC_REPTD \
  "At " MSG_TIME " repeated recoverable right-hand side function errors."
#define MSGCV_RHSFUNC_FIRST \
  "The right-hand side routine failed at the first call."
#define MSGCV_RTFUNC_FAILED                                              \
  "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable " \
  "manner."
#define MSGCV_CLOSE_ROOTS "Root found at and very near " MSG_TIME "."
#define MSGCV_BAD_TSTOP                                      \
  "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME \
  " in the direction of integration."
#define MSGCV_INACTIVE_ROOTS                                           \
  "At the end of the first step, there are still some root functions " \
  "identically 0. This warning will not be issued again."
#define MSGCV_NLS_SETUP_FAILED \
  "At " MSG_TIME ", the nonlinear solver setup failed unrecoverably."
#define MSGCV_NLS_INPUT_NULL \
  "At " MSG_TIME ", the nonlinear solver was passed a NULL input."
#define MSGCV_NLS_FAIL \
  "At " MSG_TIME ", the nonlinear solver failed in an unrecoverable manner."

/* CVode Projection Error Messages */

#define MSG_CV_MEM_NULL "cvode_mem = NULL illegal."
#define MSG_CV_MEM_FAIL "A memory request failed."

#define MSG_CV_PROJ_MEM_NULL "proj_mem = NULL illegal."
#define MSG_CV_PROJFUNC_FAIL                                              \
  "At " MSG_TIME " the projection function failed with an unrecoverable " \
  "error."
#define MSG_CV_REPTD_PROJFUNC_ERR \
  "At " MSG_TIME " the projection function had repeated recoverable errors."

#define MSGCV_NO_TOLQ \
  "No integration tolerances for quadrature variables have been specified."
#define MSGCV_BAD_EWTQ "Initial ewtQ has component(s) equal to zero (illegal)."
#define MSGCV_EWTQ_NOW_BAD \
  "At " MSG_TIME ", a component of ewtQ has become <= 0."
#define MSGCV_QRHSFUNC_FAILED                                             \
  "At " MSG_TIME ", the quadrature right-hand side routine failed in an " \
  "unrecoverable manner."
#define MSGCV_QRHSFUNC_UNREC                                                 \
  "At " MSG_TIME ", the quadrature right-hand side failed in a recoverable " \
  "manner, but no recovery is possible."
#define MSGCV_QRHSFUNC_REPTD \
  "At " MSG_TIME             \
  " repeated recoverable quadrature right-hand side function errors."
#define MSGCV_QRHSFUNC_FIRST \
  "The quadrature right-hand side routine failed at the first call."

#define MSGCV_NO_TOLS \
  "No integration tolerances for sensitivity variables have been specified."
#define MSGCV_NULL_P \
  "p = NULL when using internal DQ for sensitivity RHS illegal."
#define MSGCV_BAD_EWTS "Initial ewtS has component(s) equal to zero (illegal)."
#define MSGCV_EWTS_NOW_BAD \
  "At " MSG_TIME ", a component of ewtS has become <= 0."
#define MSGCV_SRHSFUNC_FAILED                                              \
  "At " MSG_TIME ", the sensitivity right-hand side routine failed in an " \
  "unrecoverable manner."
#define MSGCV_SRHSFUNC_UNREC                                                  \
  "At " MSG_TIME ", the sensitivity right-hand side failed in a recoverable " \
  "manner, but no recovery is possible."
#define MSGCV_SRHSFUNC_REPTD \
  "At " MSG_TIME             \
  " repeated recoverable sensitivity right-hand side function errors."
#define MSGCV_SRHSFUNC_FIRST \
  "The sensitivity right-hand side routine failed at the first call."

#define MSGCV_NULL_FQ                                                      \
  "CVODES is expected to use DQ to evaluate the RHS of quad. sensi., but " \
  "quadratures were not initialized."
#define MSGCV_NO_TOLQS                                                        \
  "No integration tolerances for quadrature sensitivity variables have been " \
  "specified."
#define MSGCV_BAD_EWTQS \
  "Initial ewtQS has component(s) equal to zero (illegal)."
#define MSGCV_EWTQS_NOW_BAD \
  "At " MSG_TIME ", a component of ewtQS has become <= 0."
#define MSGCV_QSRHSFUNC_FAILED                                           \
  "At " MSG_TIME ", the quadrature sensitivity right-hand side routine " \
  "failed in an unrecoverable manner."
#define MSGCV_QSRHSFUNC_UNREC                                                \
  "At " MSG_TIME ", the quadrature sensitivity right-hand side failed in a " \
  "recoverable manner, but no recovery is possible."
#define MSGCV_QSRHSFUNC_REPTD                                               \
  "At " MSG_TIME " repeated recoverable quadrature sensitivity right-hand " \
  "side function errors."
#define MSGCV_QSRHSFUNC_FIRST                                               \
  "The quadrature sensitivity right-hand side routine failed at the first " \
  "call."

/*
 * =================================================================
 *   A D J O I N T    E R R O R    M E S S A G E S
 * =================================================================
 */

#define MSGCV_NO_ADJ     "Illegal attempt to call before calling CVodeAdjMalloc."
#define MSGCV_BAD_STEPS  "Steps nonpositive illegal."
#define MSGCV_BAD_INTERP "Illegal value for interp."
#define MSGCV_BAD_WHICH  "Illegal value for which."
#define MSGCV_NO_BCK     "No backward problems have been defined yet."
#define MSGCV_NO_FWD     "Illegal attempt to call before calling CVodeF."
#define MSGCV_BAD_TB0                                                       \
  "The initial time tB0 for problem %d is outside the interval over which " \
  "the forward problem was solved."
#define MSGCV_BAD_SENSI                                                      \
  "At least one backward problem requires sensitivities, but they were not " \
  "stored for interpolation."
#define MSGCV_BAD_ITASKB \
  "Illegal value for itaskB. Legal values are CV_NORMAL and CV_ONE_STEP."
#define MSGCV_BAD_TBOUT                                                  \
  "The final time tBout is outside the interval over which the forward " \
  "problem was solved."
#define MSGCV_BACK_ERROR \
  "Error occurred while integrating backward problem # %d"
#define MSGCV_BAD_TINTERP "Bad t = " SUN_FORMAT_G " for interpolation."
#define MSGCV_WRONG_INTERP \
  "This function cannot be called for the specified interp type."

#ifdef __cplusplus
}
#endif

#endif
