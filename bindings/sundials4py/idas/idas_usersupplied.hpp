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

#ifndef _SUNDIALS4PY_IDAS_USERSUPPLIED_HPP
#define _SUNDIALS4PY_IDAS_USERSUPPLIED_HPP

#include <idas/idas.h>
#include <idas/idas_ls.h>

#include <sundials/sundials_core.hpp>

#include "sundials4py_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// IDAS user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to IDAS as user_data.
///////////////////////////////////////////////////////////////////////////////

struct idas_user_supplied_fn_table
{
  // user-supplied function pointers
  nb::object res, rootfn, ewtn, rwtn, resNLS;

  // idas_ls user-supplied function pointers
  nb::object lsjacfn, lsprecsetupfn, lsprecsolvefn, lsjactimessetupfn,
    lsjactimesvecfn, lsjacresfn;

  // idas quadrature user-supplied function pointers
  nb::object resQ, resQS;

  // idas FSA user-supplied function pointers
  nb::object resS;

  // idas adjoint user-supplied function pointers
  nb::object resB, resQB, resBS, resQBS;

  // idas_ls adjoint user-supplied function pointers
  nb::object lsjacfnB, lsjacfnBS, lsprecsetupfnB, lsprecsetupfnBS,
    lsprecsolvefnB, lsprecsolvefnBS, lsjactimessetupfnB, lsjactimessetupfnBS,
    lsjactimesvecfnB, lsjactimesvecfnBS;
};

inline idas_user_supplied_fn_table* get_idas_fn_table(void* ida_mem)
{
  auto mem      = static_cast<IDAMem>(ida_mem);
  auto fn_table = static_cast<idas_user_supplied_fn_table*>(mem->python);
  if (!fn_table)
  {
    throw sundials4py::null_function_table(
      "Failed to get Python function table from IDAS memory");
  }
  return fn_table;
}

///////////////////////////////////////////////////////////////////////////////
// IDAS user-supplied functions
///////////////////////////////////////////////////////////////////////////////

template<typename... Args>
inline int idas_res_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAResFn>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::res, args...);
}

using IDARootStdFn = int(sunrealtype t, N_Vector y, N_Vector yp,
                         sundials4py::Array1d gout, void* user_data);

inline int idas_rootfn_wrapper(sunrealtype t, N_Vector y, N_Vector yp,
                               sunrealtype* gout_1d, void* user_data)
{
  auto ida_mem  = static_cast<IDAMem>(user_data);
  auto fn_table = get_idas_fn_table(user_data);
  auto fn       = nb::cast<std::function<IDARootStdFn>>(fn_table->rootfn);
  auto nrtfn    = ida_mem->ida_nrtfn;

  sundials4py::Array1d gout(gout_1d, {static_cast<size_t>(nrtfn)},
                            nb::find(gout_1d));

  return fn(t, y, yp, gout, nullptr);
}

template<typename... Args>
inline int idas_ewtfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAEwtFn>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::ewtn, args...);
}

template<typename... Args>
inline int idas_nlsresfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAResFn>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::resNLS, args...);
}

template<typename... Args>
inline int idas_lsjacfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacFn>, idas_user_supplied_fn_table, IDAMem,
    4>(&idas_user_supplied_fn_table::lsjacfn, args...);
}

template<typename... Args>
inline int idas_lsprecsetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSetupFn>, idas_user_supplied_fn_table,
    IDAMem, 1>(&idas_user_supplied_fn_table::lsprecsetupfn, args...);
}

template<typename... Args>
inline int idas_lsprecsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSolveFn>, idas_user_supplied_fn_table,
    IDAMem, 1>(&idas_user_supplied_fn_table::lsprecsolvefn, args...);
}

template<typename... Args>
inline int idas_lsjactimessetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesSetupFn>, idas_user_supplied_fn_table,
    IDAMem, 1>(&idas_user_supplied_fn_table::lsjactimessetupfn, args...);
}

template<typename... Args>
inline int idas_lsjactimesvecfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesVecFn>, idas_user_supplied_fn_table,
    IDAMem, 3>(&idas_user_supplied_fn_table::lsjactimesvecfn, args...);
}

template<typename... Args>
inline int idas_lsjacresfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAResFn>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::lsjacresfn, args...);
}

template<typename... Args>
inline int idas_resQ_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAQuadRhsFn>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::resQ, args...);
}

using IDAQuadSensRhsStdFn = int(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                                std::vector<N_Vector> yyS,
                                std::vector<N_Vector> ypS, N_Vector rrQ,
                                std::vector<N_Vector> rhsvalQS, void* user_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

inline int idas_resQS_wrapper(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                              N_Vector* yyS, N_Vector* ypS, N_Vector rrQ,
                              N_Vector* rhsvalQS, void* user_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto fn_table = static_cast<idas_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<IDAQuadSensRhsStdFn>>(fn_table->resQS);

  std::vector<N_Vector> yyS_1d(yyS, yyS + Ns);
  std::vector<N_Vector> ypS_1d(ypS, ypS + Ns);
  std::vector<N_Vector> rhsvalQS_1d(rhsvalQS, rhsvalQS + Ns);

  return fn(Ns, t, yy, yp, yyS_1d, ypS_1d, rrQ, rhsvalQS_1d, nullptr, tmp1,
            tmp2, tmp3);
}

using IDASensResStdFn = int(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                            N_Vector resval, std::vector<N_Vector> yyS_1d,
                            std::vector<N_Vector> ypS_1d,
                            std::vector<N_Vector> resvalS_1d, void* user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

inline int idas_resS_wrapper(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                             N_Vector resval, N_Vector* yyS_1d, N_Vector* ypS_1d,
                             N_Vector* resvalS_1d, void* user_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto fn_table = get_idas_fn_table(user_data);
  auto fn       = nb::cast<std::function<IDASensResStdFn>>(fn_table->resS);

  std::vector<N_Vector> yyS(yyS_1d, yyS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);
  std::vector<N_Vector> resvalS(resvalS_1d, resvalS_1d + Ns);

  return fn(Ns, t, yy, yp, resval, yyS, ypS, resvalS, nullptr, tmp1, tmp2, tmp3);
}

///////////////////////////////////////////////////////////////////////////////
// IDAS Adjoint user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline idas_user_supplied_fn_table* get_idas_fn_table(void* ida_mem, int which)
{
  auto mem      = static_cast<IDAMem>(IDAGetAdjIDABmem(ida_mem, which));
  auto fn_table = static_cast<idas_user_supplied_fn_table*>(mem->python);
  if (!fn_table)
  {
    throw sundials4py::null_function_table(
      "Failed to get Python function table from IDAS memory");
  }
  return fn_table;
}

template<typename... Args>
inline int idas_resB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAResFnB>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::resB, args...);
}

template<typename... Args>
inline int idas_resQB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDAQuadRhsFnB>, idas_user_supplied_fn_table, IDAMem,
    1>(&idas_user_supplied_fn_table::resQB, args...);
}

template<typename... Args>
inline int idas_lsjacfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacFnB>, idas_user_supplied_fn_table, IDAMem,
    4>(&idas_user_supplied_fn_table::lsjacfnB, args...);
}

template<typename... Args>
inline int idas_lsprecsetupfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSetupFnB>, idas_user_supplied_fn_table,
    IDAMem, 1>(&idas_user_supplied_fn_table::lsprecsetupfnB, args...);
}

template<typename... Args>
inline int idas_lsprecsolvefnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsPrecSolveFnB>, idas_user_supplied_fn_table,
    IDAMem, 1>(&idas_user_supplied_fn_table::lsprecsolvefnB, args...);
}

template<typename... Args>
inline int idas_lsjactimessetupfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesSetupFnB>, idas_user_supplied_fn_table,
    IDAMem, 1>(&idas_user_supplied_fn_table::lsjactimessetupfnB, args...);
}

template<typename... Args>
inline int idas_lsjactimesvecfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<IDALsJacTimesVecFnB>, idas_user_supplied_fn_table,
    IDAMem, 3>(&idas_user_supplied_fn_table::lsjactimesvecfnB, args...);
}

using IDAResStdFnBS = int(sunrealtype t, N_Vector y, N_Vector yp,
                          std::vector<N_Vector> yS_1d,
                          std::vector<N_Vector> ypS_1d, N_Vector yB,
                          N_Vector ypB, N_Vector yBdot, void* user_dataB);

inline int ida_resBS_wrapper(sunrealtype t, N_Vector y, N_Vector yp,
                             N_Vector* yS_1d, N_Vector* ypS_1d, N_Vector yB,
                             N_Vector ypB, N_Vector yBdot, void* user_dataB)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = get_idas_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<IDAResStdFnBS>>(fn_table->resBS);
  auto Ns       = ida_mem->ida_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(t, y, yp, yS, ypS, yB, ypB, yBdot, nullptr);
}

using IDAQuadRhsStdFnBS = int(sunrealtype t, N_Vector y, N_Vector yp,
                              std::vector<N_Vector> yS_1d,
                              std::vector<N_Vector> ypS_1d, N_Vector yB,
                              N_Vector ypB, N_Vector rhsvalBQS, void* user_dataB);

inline int idas_resQBS_wrapper(sunrealtype t, N_Vector y, N_Vector yp,
                               N_Vector* yS_1d, N_Vector* ypS_1d, N_Vector yB,
                               N_Vector ypB, N_Vector rhsvalBQS, void* user_dataB)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = static_cast<idas_user_supplied_fn_table*>(user_dataB);
  auto fn       = nb::cast<std::function<IDAQuadRhsStdFnBS>>(fn_table->resQBS);
  auto Ns       = ida_mem->ida_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(t, y, yp, yS, ypS, yB, ypB, rhsvalBQS, nullptr);
}

using IDALsJacStdFnBS = int(sunrealtype tt, sunrealtype c_jB, N_Vector yy,
                            N_Vector yp, std::vector<N_Vector> yS_1d,
                            std::vector<N_Vector> ypS_1d, N_Vector yyB,
                            N_Vector ypB, N_Vector rrB, SUNMatrix JacB,
                            void* user_dataB, N_Vector tmp1B, N_Vector tmp2B,
                            N_Vector tmp3B);

inline int idas_lsjacfnBS_wrapper(sunrealtype tt, sunrealtype c_jB, N_Vector yy,
                                  N_Vector yp, N_Vector* yS_1d,
                                  N_Vector* ypS_1d, N_Vector yyB, N_Vector ypB,
                                  N_Vector rrB, SUNMatrix JacB, void* user_dataB,
                                  N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = get_idas_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<IDALsJacStdFnBS>>(fn_table->lsjacfnBS);
  auto Ns       = ida_mem->ida_Ns;

  std::vector<N_Vector> yS(yS_1d, yS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(tt, c_jB, yy, yp, yS, ypS, yyB, ypB, rrB, JacB, nullptr, tmp1B,
            tmp2B, tmp3B);
}

using IDALsPrecSetupStdFnBS = int(sunrealtype tt, N_Vector yy, N_Vector yp,
                                  std::vector<N_Vector> yyS_1d,
                                  std::vector<N_Vector> ypS_1d, N_Vector yyB,
                                  N_Vector ypB, N_Vector rrB, sunrealtype c_jB,
                                  void* user_dataB);

inline int idas_lsprecsetupfnBS_wrapper(sunrealtype tt, N_Vector yy, N_Vector yp,
                                        N_Vector* yyS_1d, N_Vector* ypS_1d,
                                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                        sunrealtype c_jB, void* user_dataB)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = get_idas_fn_table(user_dataB);
  auto fn =
    nb::cast<std::function<IDALsPrecSetupStdFnBS>>(fn_table->lsprecsetupfnBS);
  auto Ns = ida_mem->ida_Ns;

  std::vector<N_Vector> yyS(yyS_1d, yyS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(tt, yy, yp, yyS, ypS, yyB, ypB, rrB, c_jB, nullptr);
}

using IDALsPrecSolveStdFnBS = int(sunrealtype tt, N_Vector yy, N_Vector yp,
                                  std::vector<N_Vector> yyS_1d,
                                  std::vector<N_Vector> ypS_1d, N_Vector yyB,
                                  N_Vector ypB, N_Vector rrB, N_Vector rvecB,
                                  N_Vector zvecB, sunrealtype c_jB,
                                  sunrealtype deltaB, void* user_dataB);

inline int idas_lsprecsolvefnBS_wrapper(sunrealtype tt, N_Vector yy, N_Vector yp,
                                        N_Vector* yyS_1d, N_Vector* ypS_1d,
                                        N_Vector yyB, N_Vector ypB,
                                        N_Vector rrB, N_Vector rvecB,
                                        N_Vector zvecB, sunrealtype c_jB,
                                        sunrealtype deltaB, void* user_dataB)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = get_idas_fn_table(user_dataB);
  auto fn =
    nb::cast<std::function<IDALsPrecSolveStdFnBS>>(fn_table->lsprecsolvefnBS);
  auto Ns = ida_mem->ida_Ns;

  std::vector<N_Vector> yyS(yyS_1d, yyS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(tt, yy, yp, yyS, ypS, yyB, ypB, rrB, rvecB, zvecB, c_jB, deltaB,
            nullptr);
}

using IDALsJacTimesSetupStdFnBS = int(sunrealtype t, N_Vector yy, N_Vector yp,
                                      std::vector<N_Vector> yyS_1d,
                                      std::vector<N_Vector> ypS_1d,
                                      N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                      sunrealtype c_jB, void* user_dataB);

inline int idas_lsjactimessetupfnBS_wrapper(sunrealtype t, N_Vector yy,
                                            N_Vector yp, N_Vector* yyS_1d,
                                            N_Vector* ypS_1d, N_Vector yyB,
                                            N_Vector ypB, N_Vector rrB,
                                            sunrealtype c_jB, void* user_dataB)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = get_idas_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<IDALsJacTimesSetupStdFnBS>>(
    fn_table->lsjactimessetupfnBS);
  auto Ns = ida_mem->ida_Ns;

  std::vector<N_Vector> yyS(yyS_1d, yyS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(t, yy, yp, yyS, ypS, yyB, ypB, rrB, c_jB, nullptr);
}

using IDALsJacTimesVecStdFnBS = int(sunrealtype t, N_Vector yy, N_Vector yp,
                                    std::vector<N_Vector> yyS_1d,
                                    std::vector<N_Vector> ypS_1d, N_Vector yyB,
                                    N_Vector ypB, N_Vector rrB, N_Vector vB,
                                    N_Vector JvB, sunrealtype c_jB,
                                    void* user_dataB, N_Vector tmp1B,
                                    N_Vector tmp2B);

inline int idas_lsjactimesvecfnBS_wrapper(
  sunrealtype t, N_Vector yy, N_Vector yp, N_Vector* yyS_1d, N_Vector* ypS_1d,
  N_Vector yyB, N_Vector ypB, N_Vector rrB, N_Vector vB, N_Vector JvB,
  sunrealtype c_jB, void* user_dataB, N_Vector tmp1B, N_Vector tmp2B)
{
  auto ida_mem  = static_cast<IDAMem>(user_dataB);
  auto fn_table = get_idas_fn_table(user_dataB);
  auto fn       = nb::cast<std::function<IDALsJacTimesVecStdFnBS>>(
    fn_table->lsjactimesvecfnBS);
  auto Ns = ida_mem->ida_Ns;

  std::vector<N_Vector> yyS(yyS_1d, yyS_1d + Ns);
  std::vector<N_Vector> ypS(ypS_1d, ypS_1d + Ns);

  return fn(t, yy, yp, yyS, ypS, yyB, ypB, rrB, vB, JvB, c_jB, nullptr, tmp1B,
            tmp2B);
}

#endif
