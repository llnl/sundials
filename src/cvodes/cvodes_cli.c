/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This file provides command-line control over optional inputs
 * to CVODE.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "cvodes/cvodes.h"
#include "cvodes/cvodes_ls.h"
#include "cvodes/cvodes_proj.h"
#include "cvodes_impl.h"
#include "sundials_cli.h"

static int cvSetFromCommandLine(void* cvode_mem, const char* cvid, int argc,
                                char* argv[]);

/*---------------------------------------------------------------
  CVodeSetOptions:

  Sets CVODE options using strings.
  ---------------------------------------------------------------*/

int CVodeSetOptions(void* cvode_mem, const char* cvid, const char* file_name,
                    int argc, char* argv[])
{
  if (file_name != NULL && strlen(file_name) > 0)
  {
    int retval = CV_ILL_INPUT;
    cvProcessError(cvode_mem, retval, __LINE__, __func__, __FILE__,
                   "file-based options are not currently supported.");
    return retval;
  }

  if (argc > 0 && argv != NULL)
  {
    int retval = cvSetFromCommandLine(cvode_mem, cvid, argc, argv);
    if (retval != CV_SUCCESS) { return retval; }
  }

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to control CVODE options from the command line
 */

static int cvSetFromCommandLine(void* cvode_mem, const char* cvid, int argc,
                                char* argv[])
{
  CVodeMem cv_mem;
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  /* Set lists of command-line arguments, and the corresponding set routines */
  static const struct sunKeyIntPair int_pairs[] =
    {{"max_conv_fails", CVodeSetMaxConvFails},
     {"max_err_test_fails", CVodeSetMaxErrTestFails},
     {"max_hnil_warns", CVodeSetMaxHnilWarns},
     {"max_nonlin_iters", CVodeSetMaxNonlinIters},
     {"max_order", CVodeSetMaxOrd},
     {"stab_lim_det", CVodeSetStabLimDet},
     {"interpolate_stop_time", CVodeSetInterpolateStopTime},
     {"num_fails_eta_max_err_fail", CVodeSetNumFailsEtaMaxErrFail},
     {"quad_err_con", CVodeSetQuadErrCon},
     {"sens_err_con", CVodeSetSensErrCon},
     {"sens_max_nonlin_iters", CVodeSetSensMaxNonlinIters},
     {"quad_sens_err_con", CVodeSetQuadSensErrCon},
     {"linear_solution_scaling", CVodeSetLinearSolutionScaling},
     {"proj_err_est", CVodeSetProjErrEst},
     {"max_num_proj_fails", CVodeSetMaxNumProjFails},
     {"max_num_constraint_fails", CVodeSetMaxNumConstraintFails}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyLongPair long_pairs[] =
    {{"lsetup_frequency", CVodeSetLSetupFrequency},
     {"max_num_steps", CVodeSetMaxNumSteps},
     {"num_steps_eta_max_early_step", CVodeSetNumStepsEtaMaxEarlyStep},
     {"monitor_frequency", CVodeSetMonitorFrequency},
     {"jac_eval_frequency", CVodeSetJacEvalFrequency},
     {"proj_frequency", CVodeSetProjFrequency}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyRealPair real_pairs[] =
    {{"delta_gamma_max_lsetup", CVodeSetDeltaGammaMaxLSetup},
     {"init_step", CVodeSetInitStep},
     {"max_step", CVodeSetMaxStep},
     {"min_step", CVodeSetMinStep},
     {"nonlin_conv_coef", CVodeSetNonlinConvCoef},
     {"stop_time", CVodeSetStopTime},
     {"eta_max_first_step", CVodeSetEtaMaxFirstStep},
     {"eta_max_early_step", CVodeSetEtaMaxEarlyStep},
     {"eta_max", CVodeSetEtaMax},
     {"eta_min", CVodeSetEtaMin},
     {"eta_min_err_fail", CVodeSetEtaMinErrFail},
     {"eta_max_err_fail", CVodeSetEtaMaxErrFail},
     {"eta_conv_fail", CVodeSetEtaConvFail},
     {"delta_gamma_max_bad_jac", CVodeSetDeltaGammaMaxBadJac},
     {"eps_lin", CVodeSetEpsLin},
     {"ls_norm_factor", CVodeSetLSNormFactor},
     {"eps_proj", CVodeSetEpsProj},
     {"proj_fail_eta", CVodeSetProjFailEta}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static const struct sunKeyTwoRealPair tworeal_pairs[] =
    {{"eta_fixed_step_bounds", CVodeSetEtaFixedStepBounds},
     {"scalar_tolerances", CVodeSStolerances},
     {"quad_scalar_tolerances", CVodeQuadSStolerances}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  static const struct sunKeyTwoIntPair twoint_pairs[] =
    {{"max_order_b", CVodeSetMaxOrdB},
     {"stab_lim_det_b", CVodeSetStabLimDetB},
     {"quad_err_con_b", CVodeSetQuadErrConB},
     {"linear_solution_scaling_b", CVodeSetLinearSolutionScalingB}};
  static const int num_twoint_keys = sizeof(twoint_pairs) / sizeof(*twoint_pairs);

  static const struct sunKeyActionPair action_pairs[] =
    {{"clear_stop_time", CVodeClearStopTime},
     {"adj_no_sensi", CVodeSetAdjNoSensi},
     {"no_inactive_root_warn", CVodeSetNoInactiveRootWarn}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  static const struct sunKeyIntRealPair int_real_pairs[] =
    {{"sens_dq_method", CVodeSetSensDQMethod},
     {"init_step_b", CVodeSetInitStepB},
     {"min_step_b", CVodeSetMinStepB},
     {"max_step_b", CVodeSetMaxStepB},
     {"eps_lin_b", CVodeSetEpsLinB},
     {"ls_norm_factor_b", CVodeSetLSNormFactorB}};
  static const int num_int_real_keys = sizeof(int_real_pairs) /
                                       sizeof(*int_real_pairs);

  static const struct sunKeyIntRealRealPair int_real_real_pairs[] =
    {{"scalar_tolerances_b", CVodeSStolerancesB},
     {"quad_scalar_tolerances_b", CVodeQuadSStolerancesB}};
  static const int num_int_real_real_keys = sizeof(int_real_real_pairs) /
                                            sizeof(*int_real_real_pairs);

  static const struct sunKeyIntLongPair int_long_pairs[] = {
    {"max_num_steps_b", CVodeSetMaxNumStepsB}};
  static const int num_int_long_keys = sizeof(int_long_pairs) /
                                       sizeof(*int_long_pairs);

  /* Prefix for options to set */
  const char* default_id = "cvodes";
  size_t offset          = strlen(default_id) + 1;
  if (cvid != NULL && strlen(cvid) > 0) { offset = strlen(cvid) + 1; }
  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  if (cvid != NULL && strlen(cvid) > 0) { strcpy(prefix, cvid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  for (int idx = 1; idx < argc; idx++)
  {
    int j, retval;
    sunbooleantype arg_used = SUNFALSE;

    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* check all "int" command-line options */
    retval = sunCheckAndSetIntArgs(cvode_mem, &idx, argv, offset, int_pairs,
                                   num_int_keys, &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", int_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all long int command-line options */
    retval = sunCheckAndSetLongArgs(cvode_mem, &idx, argv, offset, long_pairs,
                                    num_long_keys, &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", long_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all real command-line options */
    retval = sunCheckAndSetRealArgs(cvode_mem, &idx, argv, offset, real_pairs,
                                    num_real_keys, &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", real_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all pair-of-int command-line options */
    retval = sunCheckAndSetTwoIntArgs(cvode_mem, &idx, argv, offset, twoint_pairs,
                                      num_twoint_keys, &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", twoint_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all pair-of-real command-line options */
    retval = sunCheckAndSetTwoRealArgs(cvode_mem, &idx, argv, offset,
                                       tworeal_pairs, num_tworeal_keys,
                                       &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", tworeal_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all action command-line options */
    retval = sunCheckAndSetActionArgs(cvode_mem, &idx, argv, offset, action_pairs,
                                      num_action_keys, &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", action_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all int+real command-line options */
    retval = sunCheckAndSetIntRealArgs(cvode_mem, &idx, argv, offset,
                                       int_real_pairs, num_int_real_keys,
                                       &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", int_real_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all int+long command-line options */
    retval = sunCheckAndSetIntLongArgs(cvode_mem, &idx, argv, offset,
                                       int_long_pairs, num_int_long_keys,
                                       &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", int_long_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all int+real+real command-line options */
    retval = sunCheckAndSetIntRealRealArgs(cvode_mem, &idx, argv, offset,
                                           int_real_real_pairs,
                                           num_int_real_real_keys, &arg_used, &j);
    if (retval != CV_SUCCESS)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "error setting key: %s", int_real_real_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* warn for uninterpreted cvid.X arguments */
    cvProcessError(cv_mem, CV_WARNING, __LINE__, __func__, __FILE__,
                   "WARNING: key %s was not handled\n", argv[idx]);
  }
  free(prefix);
  return (CV_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
