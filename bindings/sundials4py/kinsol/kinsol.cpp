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

#include "sundials4py.hpp"

#include <sundials/sundials_core.hpp>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol.hpp>
#include <kinsol/kinsol_ls.h>

#include "kinsol/kinsol_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

#include "kinsol_usersupplied.hpp"

#define BIND_KINSOL_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)       \
  m.def(                                                                \
    #NAME,                                                              \
    [](void* kin_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                   \
      auto fntable    = get_kinsol_fn_table(kin_mem);                   \
      fntable->MEMBER = nb::cast(fn);                                   \
      if (fn) { return NAME(kin_mem, WRAPPER); }                        \
      else { return NAME(kin_mem, nullptr); }                           \
    },                                                                  \
    __VA_ARGS__)

#define BIND_KINSOL_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2, \
                              MEMBER2, WRAPPER2, ...)                      \
  m.def(                                                                   \
    #NAME,                                                                 \
    [](void* kin_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,  \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                 \
    {                                                                      \
      auto fntable     = get_kinsol_fn_table(kin_mem);                     \
      fntable->MEMBER1 = nb::cast(fn1);                                    \
      fntable->MEMBER2 = nb::cast(fn2);                                    \
      if (fn1 && fn2) { return NAME(kin_mem, WRAPPER1, WRAPPER2); }        \
      else if (fn1) { return NAME(kin_mem, WRAPPER1, nullptr); }           \
      else if (fn2) { return NAME(kin_mem, nullptr, WRAPPER2); }           \
      else { return NAME(kin_mem, nullptr, nullptr); }                     \
    },                                                                     \
    __VA_ARGS__)

namespace sundials4py {

void bind_kinsol(nb::module_& m)
{
#include "kinsol_generated.hpp"

  nb::class_<KINView>(m, "KINView")
    .def("get", nb::overload_cast<>(&KINView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "KINSetOptions",
    [](void* kin_mem, const std::string& kinid, const std::string& file_name,
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

      return KINSetOptions(kin_mem, kinid.empty() ? nullptr : kinid.c_str(),
                           file_name.empty() ? nullptr : file_name.c_str(),
                           argc, argv.data());
    },
    nb::arg("kin_mem"), nb::arg("kinid"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  m.def(
    "KINCreate",
    [](SUNContext sunctx)
    { return std::make_shared<KINView>(KINCreate(sunctx)); },
    nb::arg("sunctx"), nb::keep_alive<0, 1>());

  m.def(
    "KINInit",
    [](void* kin_mem, std::function<std::remove_pointer_t<KINSysFn>> sysfn,
       N_Vector tmpl)
    {
      if (!sysfn) { throw sundials4py::illegal_value("sysfn was null"); }
      int kin_status = KINInit(kin_mem, kinsol_sysfn_wrapper, tmpl);
      if (kin_status != KIN_SUCCESS) { return kin_status; }

      auto fn_table      = new kinsol_user_supplied_fn_table;
      auto kinsol_mem    = static_cast<KINMem>(kin_mem);
      kinsol_mem->python = fn_table;
      kin_status         = KINSetUserData(kin_mem, kin_mem);
      if (kin_status != KIN_SUCCESS)
      {
        free(fn_table);
        return kin_status;
      }

      fn_table->sysfn = nb::cast(sysfn);
      return kin_status;
    },
    nb::arg("kin_mem"), nb::arg("sysfn"), nb::arg("tmpl"));

  BIND_KINSOL_CALLBACK(KINSetSysFunc, KINSysFn, sysfn, kinsol_sysfn_wrapper,
                       nb::arg("kin_mem"), nb::arg("sysfn"));

  BIND_KINSOL_CALLBACK(KINSetDampingFn, KINDampingStdFn, dampingfn,
                       kinsol_dampingfn_wrapper, nb::arg("kin_mem"),
                       nb::arg("damping_fn").none());

  BIND_KINSOL_CALLBACK(KINSetDepthFn, KINDepthStdFn, depthfn,
                       kinsol_depthfn_wrapper, nb::arg("kin_mem"),
                       nb::arg("depth_fn").none());

  BIND_KINSOL_CALLBACK2(KINSetPreconditioner, KINLsPrecSetupFn, lsprecsetupfn,
                        kinsol_lsprecsetupfn_wrapper, KINLsPrecSolveFn,
                        lsprecsolvefn, kinsol_lsprecsolvefn_wrapper,
                        nb::arg("kin_mem"), nb::arg("psetup").none(),
                        nb::arg("psolve").none());

  BIND_KINSOL_CALLBACK(KINSetJacFn, KINSysFn, lsjacfn, kinsol_lsjacfn_wrapper,
                       nb::arg("kin_mem"), nb::arg("jac").none());

  BIND_KINSOL_CALLBACK(KINSetJacTimesVecFn, KINLsJacTimesVecStdFn,
                       lsjactimesvecfn, kinsol_lsjactimesvecfn_wrapper,
                       nb::arg("kin_mem"), nb::arg("jtimes").none());

  BIND_KINSOL_CALLBACK(KINSetJacTimesVecSysFn, KINSysFn, lsjtvsysfn,
                       kinsol_lsjtvsysfn_wrapper, nb::arg("kin_mem"),
                       nb::arg("jtvSysFn").none());
}

} // namespace sundials4py

// The destroy functions gets called in our C code by the integrator destructor
extern "C" void kinsol_user_supplied_fn_table_destroy(void* ptr)
{
  delete static_cast<kinsol_user_supplied_fn_table*>(ptr);
}
