/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This file provides command-line control over optional inputs
 * to KINSOL.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_ls.h"
#include "kinsol_impl.h"
#include "sundials_cli.h"

static int kinSetFromCommandLine(void* kinmem, const char* kinid, int argc,
                                 char* argv[]);

/*---------------------------------------------------------------
  KINSetOptions:

  Sets KINSOL options using strings.
  ---------------------------------------------------------------*/

int KINSetOptions(void* kinmem, const char* kinid, const char* file_name,
                  int argc, char* argv[])
{
  if (file_name != NULL && strlen(file_name) > 0)
  {
    int retval = KIN_ILL_INPUT;
    KINProcessError(kinmem, retval, __LINE__, __func__, __FILE__,
                    "file-based options are not currently supported.");
    return retval;
  }

  if (argc > 0 && argv != NULL)
  {
    int retval = kinSetFromCommandLine(kinmem, kinid, argc, argv);
    if (retval != KIN_SUCCESS) { return retval; }
  }

  return KIN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to control KINSOL options from the command line
 */

static int kinSetFromCommandLine(void* kinmem, const char* kinid, int argc,
                                 char* argv[])
{
  KINMem kin_mem;
  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }
  kin_mem = (KINMem)kinmem;

  /* Set lists of command-line arguments, and the corresponding set routines */
  static const struct sunKeyIntPair int_pairs[] =
    {{"orth_aa", KINSetOrthAA},
     {"return_newest", KINSetReturnNewest},
     {"no_init_setup", KINSetNoInitSetup},
     {"no_res_mon", KINSetNoResMon},
     {"eta_form", KINSetEtaForm},
     {"no_min_eps", KINSetNoMinEps}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyLongPair long_pairs[] =
    {{"m_aa", KINSetMAA},
     {"delay_aa", KINSetDelayAA},
     {"num_max_iters", KINSetNumMaxIters},
     {"max_setup_calls", KINSetMaxSetupCalls},
     {"max_sub_setup_calls", KINSetMaxSubSetupCalls},
     {"max_beta_fails", KINSetMaxBetaFails}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyRealPair real_pairs[] =
    {{"damping", KINSetDamping},
     {"damping_aa", KINSetDampingAA},
     {"eta_const_value", KINSetEtaConstValue},
     {"res_mon_const_value", KINSetResMonConstValue},
     {"max_newton_step", KINSetMaxNewtonStep},
     {"rel_err_func", KINSetRelErrFunc},
     {"func_norm_tol", KINSetFuncNormTol},
     {"scaled_step_tol", KINSetScaledStepTol}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  static const struct sunKeyTwoRealPair tworeal_pairs[] =
    {{"eta_params", KINSetEtaParams}, {"res_mon_params", KINSetResMonParams}};
  static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                      sizeof(*tworeal_pairs);

  /* Prefix for options to set */
  const char* default_id = "kinsol";
  size_t offset          = strlen(default_id) + 1;
  if (kinid != NULL && strlen(kinid) > 0) { offset = strlen(kinid) + 1; }
  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  if (kinid != NULL && strlen(kinid) > 0) { strcpy(prefix, kinid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  for (int idx = 1; idx < argc; idx++)
  {
    int j, retval;
    sunbooleantype arg_used = SUNFALSE;

    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* check all "int" command-line options */
    retval = sunCheckAndSetIntArgs(kinmem, &idx, argv, offset, int_pairs,
                                   num_int_keys, &arg_used, &j);
    if (retval != KIN_SUCCESS)
    {
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s", int_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all long int command-line options */
    retval = sunCheckAndSetLongArgs(kinmem, &idx, argv, offset, long_pairs,
                                    num_long_keys, &arg_used, &j);
    if (retval != KIN_SUCCESS)
    {
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s", long_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all real command-line options */
    retval = sunCheckAndSetRealArgs(kinmem, &idx, argv, offset, real_pairs,
                                    num_real_keys, &arg_used, &j);
    if (retval != KIN_SUCCESS)
    {
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s", real_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* check all pair-of-real command-line options */
    retval = sunCheckAndSetTwoRealArgs(kinmem, &idx, argv, offset, tworeal_pairs,
                                       num_tworeal_keys, &arg_used, &j);
    if (retval != KIN_SUCCESS)
    {
      KINProcessError(kin_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key: %s", tworeal_pairs[j].key);
      free(prefix);
      return retval;
    }
    if (arg_used) continue;

    /* warn for uninterpreted kinid.X arguments */
    KINProcessError(kin_mem, KIN_WARNING, __LINE__, __func__, __FILE__,
                    "WARNING: key %s was not handled\n", argv[idx]);
  }
  free(prefix);
  return (KIN_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
