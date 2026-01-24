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

#ifndef _SUNDIALS4PY_CVODE_USERSUPPLIED_HPP
#define _SUNDIALS4PY_CVODE_USERSUPPLIED_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_ls.h>

#include <sundials/sundials_core.hpp>

#include "sundials/sundials_types.h"
#include "sundials4py_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// CVODE user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to CVODE as user_data.
///////////////////////////////////////////////////////////////////////////////

struct cvode_user_supplied_fn_table
{
  // user-supplied function pointers
  nb::object f, rootfn, ewtn, rwtn, fNLS, projfn;

  // cvode_ls user-supplied function pointers
  nb::object lsjacfn, lsprecsetupfn, lsprecsolvefn, lsjactimessetupfn,
    lsjactimesvecfn, lslinsysfn, lsjacrhsfn;

  // cvode quadrature user-supplied function pointers
  nb::object fQ, fQS;

  // cvode FSA user-supplied function pointers
  nb::object fS, fS1;

  // cvode adjoint user-supplied function pointers
  nb::object fB, fBS, fQB, fQBS;

  // cvode_ls adjoint user-supplied function pointers
  nb::object lsjacfnB, lsjacfnBS, lsprecsetupfnB, lsprecsetupfnBS,
    lsprecsolvefnB, lsprecsolvefnBS, lsjactimessetupfnB, lsjactimessetupfnBS,
    lsjactimesvecfnB, lsjactimesvecfnBS, lslinsysfnB, lslinsysfnBS;
};

// Helper to extract CVodeMem and function table
inline cvode_user_supplied_fn_table* get_cvode_fn_table(void* cv_mem)
{
  auto mem      = static_cast<CVodeMem>(cv_mem);
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(mem->python);
  if (!fn_table)
  {
    throw sundials4py::null_function_table(
      "Failed to get Python function table from CVODE memory");
  }
  return fn_table;
}

///////////////////////////////////////////////////////////////////////////////
// CVODE user-supplied functions
///////////////////////////////////////////////////////////////////////////////

template<typename... Args>
inline int cvode_f_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::f, args...);
}

using CVRootStdFn = int(sunrealtype t, N_Vector y, sundials4py::Array1d gout,
                        void* user_data);

inline int cvode_rootfn_wrapper(sunrealtype t, N_Vector y, sunrealtype* gout_1d,
                                void* user_data)
{
  auto cv_mem   = static_cast<CVodeMem>(user_data);
  auto fn_table = get_cvode_fn_table(user_data);
  auto fn       = nb::cast<std::function<CVRootStdFn>>(fn_table->rootfn);

  sundials4py::Array1d gout(gout_1d,
                            {static_cast<unsigned long>(cv_mem->cv_nrtfn)},
                            nb::find(gout_1d));

  return fn(t, y, gout, nullptr);
}

template<typename... Args>
inline int cvode_ewtfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVEwtFn>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::ewtn, args...);
}

template<typename... Args>
inline int cvode_nlsrhsfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::fNLS, args...);
}

template<typename... Args>
inline int cvode_lsjacfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFn>, cvode_user_supplied_fn_table, CVodeMem,
    4>(&cvode_user_supplied_fn_table::lsjacfn, args...);
}

using CVLsPrecSetupStdFn = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok, sunrealtype gamma,
  void* user_data);

inline int cvode_lsprecsetupfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                       sunbooleantype jok,
                                       sunbooleantype* jcurPtr,
                                       sunrealtype gamma, void* user_data)
{
  auto fn_table = get_cvode_fn_table(user_data);
  auto fn = nb::cast<std::function<CVLsPrecSetupStdFn>>(fn_table->lsprecsetupfn);

  auto result = fn(t, y, fy, jok, gamma, nullptr);

  *jcurPtr = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int cvode_lsprecsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFn>, cvode_user_supplied_fn_table,
    CVodeMem, 1>(&cvode_user_supplied_fn_table::lsprecsolvefn, args...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFn>, cvode_user_supplied_fn_table,
    CVodeMem, 1>(&cvode_user_supplied_fn_table::lsjactimessetupfn, args...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFn>, cvode_user_supplied_fn_table,
    CVodeMem, 2>(&cvode_user_supplied_fn_table::lsjactimesvecfn, args...);
}

using CVLsLinSysStdFn = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix M, sunbooleantype jok,
  sunrealtype gamma, void* user_data, N_Vector tmp1, N_Vector tmp2,
  N_Vector tmp3);

inline int cvode_lslinsysfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy,
                                    SUNMatrix M, sunbooleantype jok,
                                    sunbooleantype* jcur, sunrealtype gamma,
                                    void* user_data, N_Vector tmp1,
                                    N_Vector tmp2, N_Vector tmp3)
{
  auto fn_table = get_cvode_fn_table(user_data);
  auto fn = nb::cast<std::function<CVLsLinSysStdFn>>(fn_table->lslinsysfn);

  auto result = fn(t, y, fy, M, jok, gamma, nullptr, tmp1, tmp2, tmp3);

  *jcur = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int cvode_lsjacrhsfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::lsjacrhsfn, args...);
}

template<typename... Args>
inline int cvode_projfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVProjFn>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::projfn, args...);
}

template<typename... Args>
inline int cvode_fQ_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadRhsFn>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::fQ, args...);
}

using CVQuadSensRhsStdFn = int(int Ns, sunrealtype t, N_Vector y,
                               std::vector<N_Vector> yS_1d, N_Vector yQdot,
                               std::vector<N_Vector> yQSdot_1d, void* user_data,
                               N_Vector tmp, N_Vector tmpQ);

inline int cvode_fQS_wrapper(int Ns, sunrealtype t, N_Vector y, N_Vector* yS_1d,
                             N_Vector yQdot, N_Vector* yQSdot_1d,
                             void* user_data, N_Vector tmp, N_Vector tmpQ)
{
  auto fn_table = get_cvode_fn_table(user_data);
  auto fn       = nb::cast<std::function<CVQuadSensRhsStdFn>>(fn_table->fQS);

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);
  std::vector<N_Vector> yQSdot(yQSdot_1d, yQSdot_1d + Ns);

  return fn(Ns, t, y, yS, yQdot, yQSdot, nullptr, tmp, tmpQ);
}

using CVSensRhsStdFn = int(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                           std::vector<N_Vector> yS, std::vector<N_Vector> ySdot,
                           void* user_data, N_Vector tmp1, N_Vector tmp2);

inline int cvode_fS_wrapper(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                            N_Vector* yS, N_Vector* ySdot, void* user_data,
                            N_Vector tmp1, N_Vector tmp2)
{
  auto fn_table = get_cvode_fn_table(user_data);
  auto fn       = nb::cast<std::function<CVSensRhsStdFn>>(fn_table->fS);

  std::vector<N_Vector> yS_1d(yS, yS + Ns);
  std::vector<N_Vector> ySdot_1d(ySdot, ySdot + Ns);

  return fn(Ns, t, y, ydot, yS_1d, ySdot_1d, nullptr, tmp1, tmp2);
}

template<typename... Args>
inline int cvode_fS1_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVSensRhs1Fn>, cvode_user_supplied_fn_table, CVodeMem,
    3>(&cvode_user_supplied_fn_table::fS1, args...);
}

///////////////////////////////////////////////////////////////////////////////
// CVODE Adjoint user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline cvode_user_supplied_fn_table* get_cvode_fn_table(void* cv_mem, int which)
{
  auto cvb_mem  = static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(cv_mem, which));
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(cvb_mem->python);
  if (!fn_table)
  {
    throw sundials4py::null_function_table(
      "Failed to get Python adjoint function table from CVODE memory");
  }
  return fn_table;
}

template<typename... Args>
inline int cvode_fB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFnB>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::fB, args...);
}

template<typename... Args>
inline int cvode_fQB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadRhsFnB>, cvode_user_supplied_fn_table, CVodeMem,
    1>(&cvode_user_supplied_fn_table::fQB, args...);
}

template<typename... Args>
inline int cvode_lsjacfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFnB>, cvode_user_supplied_fn_table, CVodeMem,
    4>(&cvode_user_supplied_fn_table::lsjacfnB, args...);
}

using CVLsPrecSetupStdFnB = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB, sunbooleantype jokB,
  sunrealtype gammaB, void* user_dataB);

inline int cvode_lsprecsetupfnB_wrapper(sunrealtype t, N_Vector y, N_Vector yB,
                                        N_Vector fyB, sunbooleantype jokB,
                                        sunbooleantype* jcurPtrB,
                                        sunrealtype gammaB, void* user_dataB)
{
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn = nb::cast<std::function<CVLsPrecSetupStdFnB>>(fn_table->lsprecsetupfnB);

  auto result = fn(t, y, yB, fyB, jokB, gammaB, nullptr);

  *jcurPtrB = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int cvode_lsprecsolvefnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFnB>, cvode_user_supplied_fn_table,
    CVodeMem, 1>(&cvode_user_supplied_fn_table::lsprecsolvefnB, args...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFnB>, cvode_user_supplied_fn_table,
    CVodeMem, 1>(&cvode_user_supplied_fn_table::lsjactimessetupfnB, args...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFnB>, cvode_user_supplied_fn_table,
    CVodeMem, 2>(&cvode_user_supplied_fn_table::lsjactimesvecfnB, args...);
}

using CVLsLinSysStdFnB = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix AB,
  sunbooleantype jokB, sunrealtype gammaB, void* user_dataB, N_Vector tmp1B,
  N_Vector tmp2B, N_Vector tmp3B);

inline int cvode_lslinsysfnB_wrapper(sunrealtype t, N_Vector y, N_Vector yB,
                                     N_Vector fyB, SUNMatrix AB,
                                     sunbooleantype jokB, sunbooleantype* jcurB,
                                     sunrealtype gammaB, void* user_dataB,
                                     N_Vector tmp1B, N_Vector tmp2B,
                                     N_Vector tmp3B)
{
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn = nb::cast<std::function<CVLsLinSysStdFnB>>(fn_table->lslinsysfnB);

  auto result = fn(t, y, yB, fyB, AB, jokB, gammaB, nullptr, tmp1B, tmp2B, tmp3B);

  *jcurB = std::get<1>(result);

  return std::get<0>(result);
}

using CVRhsStdFnBS = int(sunrealtype t, N_Vector y, std::vector<N_Vector> yS_1d,
                         N_Vector yB, N_Vector yBdot, void* user_dataB);

inline int cvode_fBS_wrapper(sunrealtype t, N_Vector y, N_Vector* yS_1d,
                             N_Vector yB, N_Vector yBdot, void* user_dataB)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<CVRhsStdFnBS>>(fn_table->fBS);
  auto Ns       = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  return fn(t, y, yS, yB, yBdot, nullptr);
}

using CVQuadRhsStdFnBS = int(sunrealtype t, N_Vector y,
                             std::vector<N_Vector> yS_1d, N_Vector yB,
                             N_Vector qBdot, void* user_dataB);

inline int cvode_fQBS_wrapper(sunrealtype t, N_Vector y, N_Vector* yS_1d,
                              N_Vector yB, N_Vector qBdot, void* user_dataB)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<CVQuadRhsStdFnBS>>(fn_table->fQBS);
  auto Ns       = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  return fn(t, y, yS, yB, qBdot, nullptr);
}

using CVLsJacStdFnBS = int(sunrealtype t, N_Vector y,
                           std::vector<N_Vector> yS_1d, N_Vector yB,
                           N_Vector fyB, SUNMatrix JB, void* user_dataB,
                           N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

inline int cvode_lsjacfnBS_wrapper(sunrealtype t, N_Vector y, N_Vector* yS_1d,
                                   N_Vector yB, N_Vector fyB, SUNMatrix JB,
                                   void* user_dataB, N_Vector tmp1B,
                                   N_Vector tmp2B, N_Vector tmp3B)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<CVLsJacStdFnBS>>(fn_table->lsjacfnBS);
  auto Ns       = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  return fn(t, y, yS, yB, fyB, JB, nullptr, tmp1B, tmp2B, tmp3B);
}

using CVLsPrecSetupStdFnBS = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, std::vector<N_Vector> yS_1d, N_Vector yB,
  N_Vector fyB, sunbooleantype jokB, sunrealtype gammaB, void* user_dataB);

inline int cvode_lsprecsetupfnBS_wrapper(sunrealtype t, N_Vector y,
                                         N_Vector* yS_1d, N_Vector yB,
                                         N_Vector fyB, sunbooleantype jokB,
                                         sunbooleantype* jcurPtrB,
                                         sunrealtype gammaB, void* user_dataB)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn =
    nb::cast<std::function<CVLsPrecSetupStdFnBS>>(fn_table->lsprecsetupfnBS);
  auto Ns = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  auto result = fn(t, y, yS, yB, fyB, jokB, gammaB, nullptr);

  *jcurPtrB = std::get<1>(result);

  return std::get<0>(result);
}

using CVLsPrecSolveStdFnBS = int(sunrealtype t, N_Vector y,
                                 std::vector<N_Vector> yS_1d, N_Vector yB,
                                 N_Vector fyB, N_Vector rB, N_Vector zB,
                                 sunrealtype gammaB, sunrealtype deltaB,
                                 int lrB, void* user_dataB);

inline int cvode_lsprecsolvefnBS_wrapper(sunrealtype t, N_Vector y,
                                         N_Vector* yS_1d, N_Vector yB,
                                         N_Vector fyB, N_Vector rB, N_Vector zB,
                                         sunrealtype gammaB, sunrealtype deltaB,
                                         int lrB, void* user_dataB)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn =
    nb::cast<std::function<CVLsPrecSolveStdFnBS>>(fn_table->lsprecsolvefnBS);
  auto Ns = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  return fn(t, y, yS, yB, fyB, rB, zB, gammaB, deltaB, lrB, nullptr);
}

using CVLsJacTimesSetupStdFnBS = int(sunrealtype t, N_Vector y,
                                     std::vector<N_Vector> yS_1d, N_Vector yB,
                                     N_Vector fyB, void* user_dataB);

inline int cvode_lsjactimessetupfnBS_wrapper(sunrealtype t, N_Vector y,
                                             N_Vector* yS_1d, N_Vector yB,
                                             N_Vector fyB, void* user_dataB)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<CVLsJacTimesSetupStdFnBS>>(
    fn_table->lsjactimessetupfnBS);
  auto Ns = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  return fn(t, y, yS, yB, fyB, nullptr);
}

using CVLsJacTimesVecStdFnBS = int(N_Vector vB, N_Vector JvB, sunrealtype t,
                                   N_Vector y, std::vector<N_Vector> yS_1d,
                                   N_Vector yB, N_Vector fyB, void* user_dataB,
                                   N_Vector tmpB);

inline int cvode_lsjactimesvecfnBS_wrapper(N_Vector vB, N_Vector JvB,
                                           sunrealtype t, N_Vector y,
                                           N_Vector* yS_1d, N_Vector yB,
                                           N_Vector fyB, void* user_dataB,
                                           N_Vector tmpB)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn =
    nb::cast<std::function<CVLsJacTimesVecStdFnBS>>(fn_table->lsjactimesvecfnBS);
  auto Ns = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  return fn(vB, JvB, t, y, yS, yB, fyB, nullptr, tmpB);
}

using CVLsLinSysStdFnBS = std::tuple<int, sunbooleantype>(
  sunrealtype t, N_Vector y, std::vector<N_Vector> yS_1d, N_Vector yB,
  N_Vector fyB, SUNMatrix AB, sunbooleantype jokB, sunrealtype gammaB,
  void* user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

inline int cvode_lslinsysfnBS_wrapper(sunrealtype t, N_Vector y, N_Vector* yS_1d,
                                      N_Vector yB, N_Vector fyB, SUNMatrix AB,
                                      sunbooleantype jokB,
                                      sunbooleantype* jcurB, sunrealtype gammaB,
                                      void* user_dataB, N_Vector tmp1B,
                                      N_Vector tmp2B, N_Vector tmp3B)
{
  auto cv_mem   = static_cast<CVodeMem>(user_dataB);
  auto fn_table = get_cvode_fn_table(user_dataB);
  auto fn = nb::cast<std::function<CVLsLinSysStdFnBS>>(fn_table->lslinsysfnBS);
  auto Ns = cv_mem->cv_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);

  auto result = fn(t, y, yS, yB, fyB, AB, jokB, gammaB, nullptr, tmp1B, tmp2B,
                   tmp3B);

  *jcurB = std::get<1>(result);

  return std::get<0>(result);
}

#endif
