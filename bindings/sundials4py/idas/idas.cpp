/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
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
 *----------------------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_core.hpp>

#include <idas/idas.h>
#include <idas/idas.hpp>
#include <idas/idas_ls.h>

#include "idas/idas_impl.h"

#include "idas_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

#define BIND_IDA_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)          \
  m.def(                                                                \
    #NAME,                                                              \
    [](void* ida_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                   \
      auto fn_table    = get_idas_fn_table(ida_mem);                    \
      fn_table->MEMBER = nb::cast(fn);                                  \
      if (fn) { return NAME(ida_mem, WRAPPER); }                        \
      else { return NAME(ida_mem, nullptr); }                           \
    },                                                                  \
    __VA_ARGS__)

#define BIND_IDA_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,   \
                           MEMBER2, WRAPPER2, ...)                        \
  m.def(                                                                  \
    #NAME,                                                                \
    [](void* ida_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1, \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                \
    {                                                                     \
      auto fn_table     = get_idas_fn_table(ida_mem);                     \
      fn_table->MEMBER1 = nb::cast(fn1);                                  \
      fn_table->MEMBER2 = nb::cast(fn2);                                  \
      if (fn1 && fn2) { return NAME(ida_mem, WRAPPER1, WRAPPER2); }       \
      else if (fn1) { return NAME(ida_mem, WRAPPER1, nullptr); }          \
      else if (fn2) { return NAME(ida_mem, nullptr, WRAPPER2); }          \
      else { return NAME(ida_mem, nullptr, nullptr); }                    \
    },                                                                    \
    __VA_ARGS__)

#define BIND_IDAB_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)                    \
  m.def(                                                                           \
    #NAME,                                                                         \
    [](void* ida_mem, int which, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                              \
      void* user_data  = nullptr;                                                  \
      auto fn_table    = get_idas_fn_table(ida_mem, which);                        \
      fn_table->MEMBER = nb::cast(fn);                                             \
      if (fn) { return NAME(ida_mem, which, WRAPPER); }                            \
      else { return NAME(ida_mem, which, nullptr); }                               \
    },                                                                             \
    __VA_ARGS__)

#define BIND_IDAB_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,   \
                            MEMBER2, WRAPPER2, ...)                        \
  m.def(                                                                   \
    #NAME,                                                                 \
    [](void* ida_mem, int which,                                           \
       std::function<std::remove_pointer_t<FN_TYPE1>> fn1,                 \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                 \
    {                                                                      \
      void* user_data   = nullptr;                                         \
      auto fn_table     = get_idas_fn_table(ida_mem, which);               \
      fn_table->MEMBER1 = nb::cast(fn1);                                   \
      fn_table->MEMBER2 = nb::cast(fn2);                                   \
      if (fn1 && fn2) { return NAME(ida_mem, which, WRAPPER1, WRAPPER2); } \
      else if (fn1) { return NAME(ida_mem, which, WRAPPER1, nullptr); }    \
      else if (fn2) { return NAME(ida_mem, which, nullptr, WRAPPER2); }    \
      else { return NAME(ida_mem, which, nullptr, nullptr); }              \
    },                                                                     \
    __VA_ARGS__)

namespace sundials4py {

void bind_idas(nb::module_& m)
{
#include "idas_generated.hpp"

  nb::class_<IDAView>(m, "IDAView")
    .def("get", nb::overload_cast<>(&IDAView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "IDASetOptions",
    [](void* ida_mem, const std::string& idaid, const std::string& file_name,
       int argc, const std::vector<std::string>& args)
    {
      std::vector<char*> argv;
      argv.reserve(args.size());

      for (const auto& arg : args)
      {
        // We need a non-const char*, so we use data() and an explicit cast.
        // This is safe as long as the underlying std::string is not modified.
        argv.push_back(const_cast<char*>(arg.data()));
      }

      return IDASetOptions(ida_mem, idaid.empty() ? nullptr : idaid.c_str(),
                           file_name.empty() ? nullptr : file_name.c_str(),
                           argc, argv.data());
    },
    nb::arg("ida_mem"), nb::arg("idaid"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  m.def(
    "IDACreate",
    [](SUNContext sunctx)
    { return std::make_shared<IDAView>(IDACreate(sunctx)); },
    nb::arg("sunctx"), nb::keep_alive<0, 1>());

  m.def(
    "IDAInit",
    [](void* ida_mem, std::function<std::remove_pointer_t<IDAResFn>> res,
       sunrealtype t0, N_Vector yy0, N_Vector yp0)
    {
      if (!res) { throw sundials4py::illegal_value("res was None"); }

      int ida_status = IDAInit(ida_mem, idas_res_wrapper, t0, yy0, yp0);
      if (ida_status != IDA_SUCCESS) { return ida_status; }

      auto fn_table                        = new idas_user_supplied_fn_table;
      static_cast<IDAMem>(ida_mem)->python = fn_table;

      ida_status = IDASetUserData(ida_mem, ida_mem);
      if (ida_status != IDA_SUCCESS)
      {
        free(fn_table);
        return ida_status;
      }

      fn_table->res = nb::cast(res);

      return ida_status;
    },
    nb::arg("ida_mem"), nb::arg("res"), nb::arg("t0"), nb::arg("yy0"),
    nb::arg("yp0"));

  m.def(
    "IDARootInit",
    [](void* ida_mem, int nrtfn,
       std::function<std::remove_pointer_t<IDARootStdFn>> fn)
    {
      auto fn_table = get_idas_fn_table(ida_mem);

      if (fn)
      {
        fn_table->rootfn = nb::cast(fn);
        return IDARootInit(ida_mem, nrtfn, idas_rootfn_wrapper);
      }
      else { return IDARootInit(ida_mem, nrtfn, nullptr); }
    },
    nb::arg("ida_mem"), nb::arg("nrtfn"), nb::arg("fn").none());

  m.def(
    "IDAQuadInit",
    [](void* ida_mem, std::function<std::remove_pointer_t<IDAQuadRhsFn>> resQ,
       N_Vector yQ0)
    {
      auto fn_table = get_idas_fn_table(ida_mem);

      if (resQ)
      {
        fn_table->resQ = nb::cast(resQ);
        return IDAQuadInit(ida_mem, idas_resQ_wrapper, yQ0);
      }
      else { return IDAQuadInit(ida_mem, nullptr, yQ0); }
    },
    nb::arg("ida_mem"), nb::arg("resQ").none(), nb::arg("yQ0"));

  BIND_IDA_CALLBACK(IDAWFtolerances, IDAEwtFn, ewtn, idas_ewtfn_wrapper,
                    nb::arg("ida_mem"), nb::arg("efun").none());

  BIND_IDA_CALLBACK(IDASetNlsResFn, IDAResFn, resNLS, idas_nlsresfn_wrapper,
                    nb::arg("ida_mem"), nb::arg("res").none());

  BIND_IDA_CALLBACK(IDASetJacFn, IDALsJacFn, lsjacfn, idas_lsjacfn_wrapper,
                    nb::arg("ida_mem"), nb::arg("jac").none());

  BIND_IDA_CALLBACK2(IDASetPreconditioner, IDALsPrecSetupFn, lsprecsetupfn,
                     idas_lsprecsetupfn_wrapper, IDALsPrecSolveFn, lsprecsolvefn,
                     idas_lsprecsolvefn_wrapper, nb::arg("ida_mem"),
                     nb::arg("pset").none(), nb::arg("psolve").none());

  BIND_IDA_CALLBACK2(IDASetJacTimes, IDALsJacTimesSetupFn, lsjactimessetupfn,
                     idas_lsjactimessetupfn_wrapper, IDALsJacTimesVecFn,
                     lsjactimesvecfn, idas_lsjactimesvecfn_wrapper,
                     nb::arg("ida_mem"), nb::arg("jtsetup").none(),
                     nb::arg("jtimes").none());

  BIND_IDA_CALLBACK(IDASetJacTimesResFn, IDALsJacTimesVecFn, lsjacresfn,
                    idas_lsjacresfn_wrapper, nb::arg("ida_mem"),
                    nb::arg("jtimesResFn").none());

  //
  // Sensitivity and quadrature sensitivity bindings
  //

  m.def(
    "IDAQuadSensInit",
    [](void* ida_mem, std::function<IDAQuadSensRhsStdFn> resQS,
       std::vector<N_Vector> yQS0)
    {
      auto fn_table = get_idas_fn_table(ida_mem);

      if (resQS)
      {
        fn_table->resQS = nb::cast(resQS);
        return IDAQuadSensInit(ida_mem, idas_resQS_wrapper, yQS0.data());
      }
      else { return IDAQuadSensInit(ida_mem, nullptr, yQS0.data()); }
    },
    nb::arg("ida_mem"), nb::arg("resQS").none(), nb::arg("yQS0"));

  m.def(
    "IDASensInit",
    [](void* ida_mem, int Ns, int ism, std::function<IDASensResStdFn> resS,
       std::vector<N_Vector> yS0, std::vector<N_Vector> ypS0)
    {
      auto fn_table = get_idas_fn_table(ida_mem);
      if (resS)
      {
        fn_table->resS = nb::cast(resS);
        return IDASensInit(ida_mem, Ns, ism, idas_resS_wrapper, yS0.data(),
                           ypS0.data());
      }
      else
      {
        return IDASensInit(ida_mem, Ns, ism, nullptr, yS0.data(), ypS0.data());
      }
    },
    nb::arg("ida_mem"), nb::arg("Ns"), nb::arg("ism"), nb::arg("resS").none(),
    nb::arg("yS0"), nb::arg("ypS0"));

  ///
  // IDAS adjoint bindings
  ///

  m.def(
    "IDAInitB",
    [](void* ida_mem, int which,
       std::function<std::remove_pointer_t<IDAResFnB>> resB, sunrealtype tB0,
       N_Vector yyB0, N_Vector ypB0)
    {
      if (!resB) { throw sundials4py::illegal_value("resB was None"); }
      int ida_status = IDAInitB(ida_mem, which, idas_resB_wrapper, tB0, yyB0,
                                ypB0);
      if (ida_status != IDA_SUCCESS) { return ida_status; }

      auto fn_table    = new idas_user_supplied_fn_table;
      auto idab_mem    = static_cast<IDAMem>(IDAGetAdjIDABmem(ida_mem, which));
      idab_mem->python = fn_table;

      ida_status = IDASetUserDataB(ida_mem, which, idab_mem);
      if (ida_status != IDA_SUCCESS)
      {
        free(fn_table);
        return ida_status;
      }

      fn_table->resB = nb::cast(resB);
      return ida_status;
    },
    nb::arg("ida_mem"), nb::arg("which"), nb::arg("resB").none(),
    nb::arg("tB0"), nb::arg("yyB0"), nb::arg("ypB0"));

  m.def(
    "IDAQuadInitB",
    [](void* ida_mem, int which,
       std::function<std::remove_pointer_t<IDAQuadRhsFnB>> resQB, N_Vector yQBO)
    {
      auto fn_table = get_idas_fn_table(ida_mem, which);
      if (resQB)
      {
        fn_table->resQB = nb::cast(resQB);
        return IDAQuadInitB(ida_mem, which, idas_resQB_wrapper, yQBO);
      }
      else { return IDAQuadInitB(ida_mem, which, nullptr, yQBO); }
    },
    nb::arg("ida_mem"), nb::arg("which"), nb::arg("resQB").none(),
    nb::arg("yQBO"));

  BIND_IDAB_CALLBACK(IDASetJacFnB, IDALsJacFnB, lsjacfnB, idas_lsjacfnB_wrapper,
                     nb::arg("ida_mem"), nb::arg("which"),
                     nb::arg("jacB").none());

  BIND_IDAB_CALLBACK2(IDASetPreconditionerB, IDALsPrecSetupFnB, lsprecsetupfnB,
                      idas_lsprecsetupfnB_wrapper, IDALsPrecSolveFnB,
                      lsprecsolvefnB, idas_lsprecsolvefnB_wrapper,
                      nb::arg("ida_mem"), nb::arg("which"),
                      nb::arg("psetB").none(), nb::arg("psolveB").none());

  BIND_IDAB_CALLBACK2(IDASetJacTimesB, IDALsJacTimesSetupFnB, lsjactimessetupfnB,
                      idas_lsjactimessetupfnB_wrapper, IDALsJacTimesVecFnB,
                      lsjactimesvecfnB, idas_lsjactimesvecfnB_wrapper,
                      nb::arg("ida_mem"), nb::arg("which"),
                      nb::arg("jsetupB").none(), nb::arg("jtimesB").none());

  m.def(
    "IDAInitBS",
    [](void* ida_mem, int which, std::function<IDAResStdFnBS> resBS,
       sunrealtype tB0, N_Vector yyB0, N_Vector ypB0)
    {
      if (!resBS) { throw sundials4py::illegal_value("resBS was null"); }
      int ida_status = IDAInitBS(ida_mem, which, ida_resBS_wrapper, tB0, yyB0,
                                 ypB0);
      if (ida_status != IDA_SUCCESS) { return ida_status; }

      auto fn_table    = new idas_user_supplied_fn_table;
      auto idab_mem    = static_cast<IDAMem>(IDAGetAdjIDABmem(ida_mem, which));
      idab_mem->python = fn_table;

      ida_status = IDASetUserDataB(ida_mem, which, idab_mem);
      if (ida_status != IDA_SUCCESS)
      {
        free(fn_table);
        return ida_status;
      }

      fn_table->resBS = nb::cast(resBS);
      return ida_status;
    },
    nb::arg("ida_mem"), nb::arg("which"), nb::arg("resBS").none(),
    nb::arg("tB0"), nb::arg("yyB0"), nb::arg("ypB0"));

  m.def(
    "IDAQuadInitBS",
    [](void* ida_mem, int which, std::function<IDAQuadRhsStdFnBS> resQBS,
       N_Vector yQBO)
    {
      auto fn_table = get_idas_fn_table(ida_mem, which);
      if (resQBS)
      {
        fn_table->resQBS = nb::cast(resQBS);
        return IDAQuadInitBS(ida_mem, which, idas_resQBS_wrapper, yQBO);
      }
      else { return IDAQuadInitBS(ida_mem, which, nullptr, yQBO); }
    },
    nb::arg("ida_mem"), nb::arg("which"), nb::arg("resQBS").none(),
    nb::arg("yQBO"));

  BIND_IDAB_CALLBACK(IDASetJacFnBS, IDALsJacStdFnBS, lsjacfnBS,
                     idas_lsjacfnBS_wrapper, nb::arg("ida_mem"),
                     nb::arg("which"), nb::arg("jacBS").none());

  BIND_IDAB_CALLBACK2(IDASetPreconditionerBS, IDALsPrecSetupStdFnBS,
                      lsprecsetupfnBS, idas_lsprecsetupfnBS_wrapper,
                      IDALsPrecSolveStdFnBS, lsprecsolvefnBS,
                      idas_lsprecsolvefnBS_wrapper, nb::arg("ida_mem"),
                      nb::arg("which"), nb::arg("psetBS").none(),
                      nb::arg("psolveBS").none());

  BIND_IDAB_CALLBACK2(IDASetJacTimesBS, IDALsJacTimesSetupStdFnBS,
                      lsjactimessetupfnBS, idas_lsjactimessetupfnBS_wrapper,
                      IDALsJacTimesVecStdFnBS, lsjactimesvecfnBS,
                      idas_lsjactimesvecfnBS_wrapper, nb::arg("ida_mem"),
                      nb::arg("which"), nb::arg("jsetupBS").none(),
                      nb::arg("jtimesBS").none());
}

} // namespace sundials4py

// The destroy functions gets called in our C code by the integrator destructor
extern "C" void idas_user_supplied_fn_table_destroy(void* ptr)
{
  delete static_cast<idas_user_supplied_fn_table*>(ptr);
}
