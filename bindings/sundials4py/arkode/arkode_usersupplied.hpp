/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
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
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_ARKODE_USERSUPPLIED_HPP
#define _SUNDIALS4PY_ARKODE_USERSUPPLIED_HPP

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_lsrkstep.h>
#include <arkode/arkode_mristep.h>

#include <sundials/sundials_core.hpp>

#include "arkode_mristep_impl.h"

#include "sundials/sundials_nvector.h"
#include "sundials4py_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// ARKODE user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to ARKODE as user_data.
///////////////////////////////////////////////////////////////////////////////

struct arkode_user_supplied_fn_table
{
  // common user-supplied function pointers
  nb::object rootfn;
  nb::object ewtn;
  nb::object rwtn;
  nb::object vecresizefn;
  nb::object postprocessstepfn;
  nb::object postprocessstagefn;
  nb::object stagepredictfn;
  nb::object relaxfn;
  nb::object relaxjacfn;
  nb::object nlsfi;

  // arkode_ls user-supplied function pointers
  nb::object lsjacfn;
  nb::object lsmassfn;
  nb::object lsprecsetupfn;
  nb::object lsprecsolvefn;
  nb::object lsjactimessetupfn;
  nb::object lsjactimesvecfn;
  nb::object lslinsysfn;
  nb::object lsmasstimessetupfn;
  nb::object lsmasstimesvecfn;
  nb::object lsmassprecsetupfn;
  nb::object lsmassprecsolvefn;
  nb::object lsjacrhsfn;

  // erkstep-specific user-supplied function pointers
  nb::object erkstep_f;
  nb::object erkstep_adjf;

  // arkstep-specific user-supplied function pointers
  nb::object arkstep_fe;
  nb::object arkstep_fi;
  nb::object arkstep_adjfe;
  nb::object arkstep_adjfi;

  // sprkstep-specific user-supplied function pointers
  nb::object sprkstep_f1;
  nb::object sprkstep_f2;

  // lsrkstep-specific user-supplied function pointers
  nb::object lsrkstep_f;
  nb::object lsrkstep_domeig;

  // mristep-specific user-supplied function pointers
  nb::object mristep_fse;
  nb::object mristep_fsi;
  nb::object mristep_preinnerfn;
  nb::object mristep_postinnerfn;
};

struct mristepinnerstepper_user_supplied_fn_table
{
  nb::object mristepinner_evolvefn;
  nb::object mristepinner_fullrhsfn;
  nb::object mristepinner_resetfn;
  nb::object mristepinner_getaccumulatederrorfn;
  nb::object mristepinner_resetaccumulatederrorfn;
  nb::object mristepinner_setrtolfn;
};

///////////////////////////////////////////////////////////////////////////////
// ARKODE user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline arkode_user_supplied_fn_table* get_arkode_fn_table(void* ark_mem)
{
  auto mem      = static_cast<ARKodeMem>(ark_mem);
  auto fn_table = static_cast<arkode_user_supplied_fn_table*>(mem->python);
  if (!fn_table)
  {
    throw sundials4py::null_function_table(
      "Failed to get Python function table from ARKODE memory");
  }
  return fn_table;
}

using ARKRootStdFn = int(sunrealtype t, N_Vector y, sundials4py::Array1d gout,
                         void* user_data);

inline int arkode_rootfn_wrapper(sunrealtype t, N_Vector y,
                                 sunrealtype* gout_1d, void* user_data)
{
  auto fn_table = get_arkode_fn_table(user_data);
  auto fn       = nb::cast<std::function<ARKRootStdFn>>(fn_table->rootfn);
  auto nrtfn    = static_cast<ARKodeMem>(user_data)->root_mem->nrtfn;

  sundials4py::Array1d gout(gout_1d, {static_cast<unsigned long>(nrtfn)},
                            nb::find(gout_1d));

  int status = fn(t, y, gout, nullptr);

  return status;
}

template<typename... Args>
inline int arkode_ewtfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKEwtFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::ewtn, args...);
}

template<typename... Args>
inline int arkode_rwtfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRwtFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::rwtn, args...);
}

template<typename... Args>
inline int arkode_vecresizefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKVecResizeFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::vecresizefn, args...);
}

template<typename... Args>
inline int arkode_postprocessstepfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKPostProcessFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::postprocessstepfn, args...);
}

template<typename... Args>
inline int arkode_postprocessstagefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKPostProcessFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::postprocessstagefn, args...);
}

template<typename... Args>
inline int arkode_stagepredictfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKStagePredictFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::stagepredictfn, args...);
}

template<typename... Args>
inline int arkode_nlsrhsfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::nlsfi, args...);
}

using ARKRelaxStdFn = std::tuple<int, sunrealtype>(N_Vector y, void* user_data);

inline int arkode_relaxfn_wrapper(N_Vector y, sunrealtype* r, void* user_data)
{
  auto fn_table = static_cast<arkode_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<ARKRelaxStdFn>>(fn_table->relaxfn);

  auto result = fn(y, nullptr);

  *r = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int arkode_relaxjacfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRelaxJacFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::relaxjacfn, args...);
}

template<typename... Args>
inline int arkode_lsjacfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsJacFn>, arkode_user_supplied_fn_table, ARKodeMem,
    4>(&arkode_user_supplied_fn_table::lsjacfn, args...);
}

template<typename... Args>
inline int arkode_lsmassfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 4>(&arkode_user_supplied_fn_table::lsmassfn, args...);
}

using ARKLsPrecSetupStdFn = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok, sunrealtype gamma,
  void* user_data);

inline int arkode_lsprecsetupfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                        sunbooleantype jok,
                                        sunbooleantype* jcurPtr,
                                        sunrealtype gamma, void* user_data)
{
  auto fn_table = static_cast<arkode_user_supplied_fn_table*>(user_data);
  auto fn = nb::cast<std::function<ARKLsPrecSetupStdFn>>(fn_table->lsprecsetupfn);

  auto result = fn(t, y, fy, jok, gamma, nullptr);

  *jcurPtr = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int arkode_lsprecsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsPrecSolveFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::lsprecsolvefn, args...);
}

template<typename... Args>
inline int arkode_lsjactimessetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsJacTimesSetupFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::lsjactimessetupfn, args...);
}

template<typename... Args>
inline int arkode_lsjactimesvecfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsJacTimesVecFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 2>(&arkode_user_supplied_fn_table::lsjactimesvecfn, args...);
}

using ARKLsLinSysStdFn = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A, SUNMatrix M,
  sunbooleantype jok, sunrealtype gamma, void* user_data, N_Vector tmp1,
  N_Vector tmp2, N_Vector tmp3);

inline int arkode_lslinsysfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                     SUNMatrix A, SUNMatrix M,
                                     sunbooleantype jok, sunbooleantype* jcurPtr,
                                     sunrealtype gamma, void* user_data,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto fn_table = static_cast<arkode_user_supplied_fn_table*>(user_data);
  auto fn = nb::cast<std::function<ARKLsLinSysStdFn>>(fn_table->lslinsysfn);

  auto result = fn(t, y, fy, A, M, jok, gamma, nullptr, tmp1, tmp2, tmp3);

  *jcurPtr = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int arkode_lsmasstimessetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassTimesSetupFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::lsmasstimessetupfn, args...);
}

template<typename... Args>
inline int arkode_lsmasstimesvecfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassTimesVecFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::lsmasstimesvecfn, args...);
}

template<typename... Args>
inline int arkode_lsmassprecsetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassPrecSetupFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::lsmassprecsetupfn, args...);
}

template<typename... Args>
inline int arkode_lsmassprecsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKLsMassPrecSolveFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::lsmassprecsolvefn, args...);
}

template<typename... Args>
inline int arkode_lsjacrhsfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::lsjacrhsfn, args...);
}

///////////////////////////////////////////////////////////////////////////////
// ERKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

template<typename... Args>
inline int erkstep_f_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::erkstep_f, args...);
}

template<typename... Args>
inline int erkstep_adjf_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNAdjRhsFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::erkstep_adjf, args...);
}

///////////////////////////////////////////////////////////////////////////////
// ARKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline int arkstep_fe_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                              void* user_data)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::arkstep_fe, t, y, ydot, user_data);
}

inline int arkstep_fi_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                              void* user_data)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::arkstep_fi, t, y, ydot, user_data);
}

inline int arkstep_adjfe_wrapper(sunrealtype t, N_Vector y, N_Vector sens,
                                 N_Vector sens_dot, void* user_data)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNAdjRhsFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::arkstep_adjfe, t, y, sens,
                  sens_dot, user_data);
}

inline int arkstep_adjfi_wrapper(sunrealtype t, N_Vector y, N_Vector sens,
                                 N_Vector sens_dot, void* user_data)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNAdjRhsFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::arkstep_adjfi, t, y, sens,
                  sens_dot, user_data);
}

///////////////////////////////////////////////////////////////////////////////
// SPRKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

template<typename... Args>
inline int sprkstep_f1_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::sprkstep_f1, args...);
}

template<typename... Args>
inline int sprkstep_f2_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::sprkstep_f2, args...);
}

///////////////////////////////////////////////////////////////////////////////
// LSRKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

template<typename... Args>
inline int lsrkstep_f_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::lsrkstep_f, args...);
}

using ARKDomEigStdFn = std::tuple<int, sunrealtype, sunrealtype>(
  sunrealtype t, N_Vector y, N_Vector fn, void* user_data, N_Vector temp1,
  N_Vector temp2, N_Vector temp3);

inline int lsrkstep_domeig_wrapper(sunrealtype t, N_Vector y, N_Vector fn,
                                   sunrealtype* lambdaR, sunrealtype* lambdaI,
                                   void* user_data, N_Vector temp1,
                                   N_Vector temp2, N_Vector temp3)
{
  auto fn_table = get_arkode_fn_table(user_data);
  auto callback =
    nb::cast<std::function<ARKDomEigStdFn>>(fn_table->lsrkstep_domeig);

  auto result = callback(t, y, fn, nullptr, temp1, temp2, temp3);

  *lambdaR = std::get<1>(result);
  *lambdaI = std::get<2>(result);

  return std::get<0>(result);
}

///////////////////////////////////////////////////////////////////////////////
// MRIStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline int mristep_fse_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                               void* user_data)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::mristep_fse, t, y, ydot, user_data);
}

inline int mristep_fsi_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                               void* user_data)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<ARKRhsFn>, arkode_user_supplied_fn_table, ARKodeMem,
    1>(&arkode_user_supplied_fn_table::mristep_fsi, t, y, ydot, user_data);
}

using MRIStepPreInnerStdFn = int(sunrealtype t, std::vector<N_Vector> f,
                                 int nvecs, void* user_data);

inline int mristep_preinnerfn_wrapper(sunrealtype t, N_Vector* f_1d, int nvecs,
                                      void* user_data)
{
  auto fn_table = static_cast<arkode_user_supplied_fn_table*>(user_data);
  auto fn =
    nb::cast<std::function<MRIStepPreInnerStdFn>>(fn_table->mristep_preinnerfn);

  std::vector<N_Vector> f(f_1d, f_1d + nvecs);

  return fn(t, f, nvecs, nullptr);
}

template<typename... Args>
inline int mristep_postinnerfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<MRIStepPostInnerFn>, arkode_user_supplied_fn_table,
    ARKodeMem, 1>(&arkode_user_supplied_fn_table::mristep_postinnerfn, args...);
}

inline int mristepinner_evolvefn_wrapper(MRIStepInnerStepper stepper,
                                         sunrealtype t0, sunrealtype tout,
                                         N_Vector y)
{
  void* user_data = nullptr;
  MRIStepInnerStepper_GetContent(stepper, &user_data);

  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<MRIStepInnerEvolveFn>, mristepinnerstepper_user_supplied_fn_table,
    MRIStepInnerStepper>(&mristepinnerstepper_user_supplied_fn_table::mristepinner_evolvefn,
                         stepper, t0, tout, y);
}

inline int mristepinner_fullrhsfn_wrapper(MRIStepInnerStepper stepper,
                                          sunrealtype t, N_Vector y, N_Vector f,
                                          int mode)
{
  void* user_data = nullptr;
  MRIStepInnerStepper_GetContent(stepper, &user_data);

  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<MRIStepInnerFullRhsFn>, mristepinnerstepper_user_supplied_fn_table,
    MRIStepInnerStepper>(&mristepinnerstepper_user_supplied_fn_table::mristepinner_fullrhsfn,
                         stepper, t, y, f, mode);
}

inline int mristepinner_resetfn_wrapper(MRIStepInnerStepper stepper,
                                        sunrealtype tR, N_Vector yR)
{
  void* user_data = nullptr;
  MRIStepInnerStepper_GetContent(stepper, &user_data);

  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<MRIStepInnerResetFn>, mristepinnerstepper_user_supplied_fn_table,
    MRIStepInnerStepper>(&mristepinnerstepper_user_supplied_fn_table::mristepinner_resetfn,
                         stepper, tR, yR);
}

using MRIStepInnerGetAccumulatedErrorStdFn =
  std::tuple<int, sunrealtype>(MRIStepInnerStepper stepper);

inline int mristepinner_getaccumulatederrorfn_wrapper(MRIStepInnerStepper stepper,
                                                      sunrealtype* accum_error)
{
  void* user_data = nullptr;
  MRIStepInnerStepper_GetContent(stepper, &user_data);

  auto fn_table =
    static_cast<mristepinnerstepper_user_supplied_fn_table*>(user_data);
  auto fn = nb::cast<std::function<MRIStepInnerGetAccumulatedErrorStdFn>>(
    fn_table->mristepinner_getaccumulatederrorfn);

  auto result = fn(stepper);

  *accum_error = std::get<1>(result);

  return std::get<0>(result);
}

inline int mristepinner_resetaccumulatederrorfn_wrapper(MRIStepInnerStepper stepper)
{
  void* user_data = nullptr;
  MRIStepInnerStepper_GetContent(stepper, &user_data);

  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<MRIStepInnerResetAccumulatedError>,
    mristepinnerstepper_user_supplied_fn_table,
    MRIStepInnerStepper>(&mristepinnerstepper_user_supplied_fn_table::mristepinner_resetaccumulatederrorfn,
                         stepper);
}

inline int mristepinner_setrtolfn_wrapper(MRIStepInnerStepper stepper,
                                          sunrealtype rtol)
{
  void* user_data = nullptr;
  MRIStepInnerStepper_GetContent(stepper, &user_data);

  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<MRIStepInnerSetRTol>, mristepinnerstepper_user_supplied_fn_table,
    MRIStepInnerStepper>(&mristepinnerstepper_user_supplied_fn_table::mristepinner_setrtolfn,
                         stepper, rtol);
}

#endif
