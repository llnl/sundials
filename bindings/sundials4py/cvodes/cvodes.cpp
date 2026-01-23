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

#include <cvodes/cvodes.h>
#include <cvodes/cvodes.hpp>
#include <cvodes/cvodes_ls.h>
#include <cvodes/cvodes_proj.h>

#include "cvodes/cvodes_impl.h"
#include "cvodes_usersupplied.hpp"

#include "sundials_adjointcheckpointscheme_impl.h"

namespace sundials4py {

using namespace sundials::experimental;

#define BIND_CVODE_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)       \
  m.def(                                                               \
    #NAME,                                                             \
    [](void* cv_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                  \
      auto fn_table    = get_cvode_fn_table(cv_mem);                   \
      fn_table->MEMBER = nb::cast(fn);                                 \
      if (fn) { return NAME(cv_mem, WRAPPER); }                        \
      else { return NAME(cv_mem, nullptr); }                           \
    },                                                                 \
    __VA_ARGS__)

#define BIND_CVODE_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2, \
                             MEMBER2, WRAPPER2, ...)                      \
  m.def(                                                                  \
    #NAME,                                                                \
    [](void* cv_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,  \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                \
    {                                                                     \
      auto fn_table     = get_cvode_fn_table(cv_mem);                     \
      fn_table->MEMBER1 = nb::cast(fn1);                                  \
      fn_table->MEMBER2 = nb::cast(fn2);                                  \
      if (fn1 && fn2) { return NAME(cv_mem, WRAPPER1, WRAPPER2); }        \
      else if (fn1) { return NAME(cv_mem, WRAPPER1, nullptr); }           \
      else if (fn2) { return NAME(cv_mem, nullptr, WRAPPER2); }           \
      else { return NAME(cv_mem, nullptr, nullptr); }                     \
    },                                                                    \
    __VA_ARGS__)

#define BIND_CVODEB_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)                 \
  m.def(                                                                          \
    #NAME,                                                                        \
    [](void* cv_mem, int which, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                             \
      void* user_data  = nullptr;                                                 \
      auto fn_table    = get_cvode_fn_table(cv_mem, which);                       \
      fn_table->MEMBER = nb::cast(fn);                                            \
      if (fn) { return NAME(cv_mem, which, WRAPPER); }                            \
      else { return NAME(cv_mem, which, nullptr); }                               \
    },                                                                            \
    __VA_ARGS__)

#define BIND_CVODEB_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2, \
                              MEMBER2, WRAPPER2, ...)                      \
  m.def(                                                                   \
    #NAME,                                                                 \
    [](void* cv_mem, int which,                                            \
       std::function<std::remove_pointer_t<FN_TYPE1>> fn1,                 \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                 \
    {                                                                      \
      void* user_data   = nullptr;                                         \
      auto fn_table     = get_cvode_fn_table(cv_mem, which);               \
      fn_table->MEMBER1 = nb::cast(fn1);                                   \
      fn_table->MEMBER2 = nb::cast(fn2);                                   \
      if (fn1 && fn2) { return NAME(cv_mem, which, WRAPPER1, WRAPPER2); }  \
      else if (fn1) { return NAME(cv_mem, which, WRAPPER1, nullptr); }     \
      else if (fn2) { return NAME(cv_mem, which, nullptr, WRAPPER2); }     \
      else { return NAME(cv_mem, which, nullptr, nullptr); }               \
    },                                                                     \
    __VA_ARGS__)

void bind_cvodes(nb::module_& m)
{
#include "cvodes_generated.hpp"

  nb::class_<CVodeView>(m, "CVodeView")
    .def("get", nb::overload_cast<>(&CVodeView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "CVodeSetOptions",
    [](void* cv_mem, const std::string& cvid, const std::string& file_name,
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

      return CVodeSetOptions(cv_mem, cvid.empty() ? nullptr : cvid.c_str(),
                             file_name.empty() ? nullptr : file_name.c_str(),
                             argc, argv.data());
    },
    nb::arg("cv_mem"), nb::arg("cvid"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  m.def(
    "CVodeCreate",
    [](int lmm, SUNContext sunctx)
    { return std::make_shared<CVodeView>(CVodeCreate(lmm, sunctx)); },
    nb::arg("lmm"), nb::arg("sunctx"), nb::keep_alive<0, 2>());

  m.def(
    "CVodeInit",
    [](void* cv_mem, std::function<std::remove_pointer_t<CVRhsFn>> rhs,
       sunrealtype t0, N_Vector y0)
    {
      if (!rhs) { throw sundials4py::illegal_value("rhs was None"); }

      int cv_status = CVodeInit(cv_mem, cvode_f_wrapper, t0, y0);
      if (cv_status != CV_SUCCESS) { return cv_status; }

      // Create the user-supplied function table to store the Python user functions
      auto fn_table = new cvode_user_supplied_fn_table;

      static_cast<CVodeMem>(cv_mem)->python = fn_table;

      // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
      cv_status = CVodeSetUserData(cv_mem, cv_mem);
      if (cv_status != CV_SUCCESS)
      {
        free(fn_table);
        return cv_status;
      }

      // Finally, set the RHS function
      fn_table->f = nb::cast(rhs);

      return cv_status;
    },
    nb::arg("cv_mem"), nb::arg("rhs"), nb::arg("t0"), nb::arg("y0"));

  m.def(
    "CVodeRootInit",
    [](void* cv_mem, int nrtfn, std::function<std::remove_pointer_t<CVRootFn>> fn)
    {
      auto fn_table = get_cvode_fn_table(cv_mem);
      if (fn)
      {
        fn_table->rootfn = nb::cast(fn);
        return CVodeRootInit(cv_mem, nrtfn, cvode_rootfn_wrapper);
      }
      else { return CVodeRootInit(cv_mem, nrtfn, nullptr); }
    },
    nb::arg("cv_mem"), nb::arg("nrtfn"), nb::arg("fn").none());

  m.def(
    "CVodeQuadInit",
    [](void* cv_mem, std::function<std::remove_pointer_t<CVQuadRhsFn>> fQ,
       N_Vector yQ0)
    {
      auto fn_table = get_cvode_fn_table(cv_mem);
      if (fQ)
      {
        fn_table->fQ = nb::cast(fQ);
        return CVodeQuadInit(cv_mem, cvode_fQ_wrapper, yQ0);
      }
      else { return CVodeQuadInit(cv_mem, nullptr, yQ0); }
    },
    nb::arg("cv_mem"), nb::arg("fQ").none(), nb::arg("yQ0"));

  BIND_CVODE_CALLBACK(CVodeWFtolerances, CVEwtFn, ewtn, cvode_ewtfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("efun").none());

  BIND_CVODE_CALLBACK(CVodeSetNlsRhsFn, CVRhsFn, fNLS, cvode_nlsrhsfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("f").none());

  BIND_CVODE_CALLBACK(CVodeSetJacFn, CVLsJacFn, lsjacfn, cvode_lsjacfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("jac").none());

  BIND_CVODE_CALLBACK2(CVodeSetPreconditioner, CVLsPrecSetupFn, lsprecsetupfn,
                       cvode_lsprecsetupfn_wrapper, CVLsPrecSolveFn,
                       lsprecsolvefn, cvode_lsprecsolvefn_wrapper,
                       nb::arg("cvode_mem"), nb::arg("pset").none(),
                       nb::arg("psolve").none());

  BIND_CVODE_CALLBACK2(CVodeSetJacTimes, CVLsJacTimesSetupFn, lsjactimessetupfn,
                       cvode_lsjactimessetupfn_wrapper, CVLsJacTimesVecFn,
                       lsjactimesvecfn, cvode_lsjactimesvecfn_wrapper,
                       nb::arg("cvode_mem"), nb::arg("jtsetup").none(),
                       nb::arg("jtimes").none());

  BIND_CVODE_CALLBACK(CVodeSetLinSysFn, CVLsLinSysFn, lslinsysfn,
                      cvode_lslinsysfn_wrapper, nb::arg("cvode_mem"),
                      nb::arg("linsys").none());

  BIND_CVODE_CALLBACK(CVodeSetJacTimesRhsFn, CVLsJacTimesVecFn, lsjacrhsfn,
                      cvode_lsjacrhsfn_wrapper, nb::arg("cvode_mem"),
                      nb::arg("jtimesRhsFn").none());

  BIND_CVODE_CALLBACK(CVodeSetProjFn, CVProjFn, projfn, cvode_projfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("pfun").none());

  m.def(
    "CVodeQuadSensInit",
    [](void* cv_mem, std::function<CVQuadSensRhsStdFn> fQS,
       std::vector<N_Vector> yQS0)
    {
      auto fn_table = get_cvode_fn_table(cv_mem);
      if (fQS)
      {
        fn_table->fQS = nb::cast(fQS);
        return CVodeQuadSensInit(cv_mem, cvode_fQS_wrapper, yQS0.data());
      }
      else { return CVodeQuadSensInit(cv_mem, nullptr, yQS0.data()); }
    },
    nb::arg("cvode_mem"), nb::arg("fQS"), nb::arg("yQS0"));

  m.def(
    "CVodeSensInit",
    [](void* cv_mem, int Ns, int ism, std::function<CVSensRhsStdFn> fS,
       std::vector<N_Vector> yS0)
    {
      auto fn_table = get_cvode_fn_table(cv_mem);

      if (fS)
      {
        fn_table->fS = nb::cast(fS);
        return CVodeSensInit(cv_mem, Ns, ism, cvode_fS_wrapper, yS0.data());
      }
      else { return CVodeSensInit(cv_mem, Ns, ism, nullptr, yS0.data()); }
    },
    nb::arg("cvode_mem"), nb::arg("Ns"), nb::arg("ism"), nb::arg("fS").none(),
    nb::arg("yS0"));

  m.def(
    "CVodeSensInit1",
    [](void* cv_mem, int Ns, int ism,
       std::function<std::remove_pointer_t<CVSensRhs1Fn>> fS1,
       std::vector<N_Vector> yS0)
    {
      auto fn_table = get_cvode_fn_table(cv_mem);

      if (fS1)
      {
        fn_table->fS1 = nb::cast(fS1);
        return CVodeSensInit1(cv_mem, Ns, ism, cvode_fS1_wrapper, yS0.data());
      }
      else { return CVodeSensInit1(cv_mem, Ns, ism, nullptr, yS0.data()); }
    },
    nb::arg("cvode_mem"), nb::arg("Ns"), nb::arg("ism"), nb::arg("fS1").none(),
    nb::arg("yS0"));

  ///
  // CVODES Adjoint Bindings
  ///

  m.def(
    "CVodeInitB",
    [](void* cv_mem, int which, std::function<std::remove_pointer_t<CVRhsFnB>> fB,
       sunrealtype tB0, N_Vector yB0)
    {
      if (!fB) { throw sundials4py::illegal_value("rhs was None"); }

      int cv_status = CVodeInitB(cv_mem, which, cvode_fB_wrapper, tB0, yB0);
      if (cv_status != CV_SUCCESS) { return cv_status; }

      auto fn_table = new cvode_user_supplied_fn_table;
      auto cvb_mem = static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(cv_mem, which));
      cvb_mem->python = fn_table;

      cv_status = CVodeSetUserDataB(cv_mem, which, cvb_mem);
      if (cv_status != CV_SUCCESS)
      {
        free(fn_table);
        return cv_status;
      }

      fn_table->fB = nb::cast(fB);
      return cv_status;
    },
    nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fB").none(),
    nb::arg("tB0"), nb::arg("yB0"));

  m.def(
    "CVodeQuadInitB",
    [](void* cv_mem, int which,
       std::function<std::remove_pointer_t<CVQuadRhsFnB>> fQB, N_Vector yQBO)
    {
      auto fn_table = get_cvode_fn_table(cv_mem, which);

      if (fQB)
      {
        fn_table->fQB = nb::cast(fQB);
        return CVodeQuadInitB(cv_mem, which, cvode_fQB_wrapper, yQBO);
      }
      else { return CVodeQuadInitB(cv_mem, which, nullptr, yQBO); }
    },
    nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fQB").none(),
    nb::arg("yQB0"));

  BIND_CVODEB_CALLBACK(CVodeSetJacFnB, CVLsJacFnB, lsjacfnB,
                       cvode_lsjacfnB_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("jacB").none());

  BIND_CVODEB_CALLBACK2(CVodeSetPreconditionerB, CVLsPrecSetupFnB, lsprecsetupfnB,
                        cvode_lsprecsetupfnB_wrapper, CVLsPrecSolveFnB,
                        lsprecsolvefnB, cvode_lsprecsolvefnB_wrapper,
                        nb::arg("cv_mem"), nb::arg("which"),
                        nb::arg("psetB").none(), nb::arg("psolveB").none());

  BIND_CVODEB_CALLBACK2(CVodeSetJacTimesB, CVLsJacTimesSetupFnB,
                        lsjactimessetupfnB, cvode_lsjactimessetupfnB_wrapper,
                        CVLsJacTimesVecFnB, lsjactimesvecfnB,
                        cvode_lsjactimesvecfnB_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("jsetupB").none(),
                        nb::arg("jtimesB").none());

  BIND_CVODEB_CALLBACK(CVodeSetLinSysFnB, CVLsLinSysFnB, lslinsysfnB,
                       cvode_lslinsysfnB_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("linsysB").none());

  m.def(
    "CVodeInitBS",
    [](void* cv_mem, int which, std::function<CVRhsStdFnBS> fBS,
       sunrealtype tB0, N_Vector yB0)
    {
      if (!fBS) { throw sundials4py::illegal_value("rhs was None"); }

      int cv_status = CVodeInitBS(cv_mem, which, cvode_fBS_wrapper, tB0, yB0);
      if (cv_status != CV_SUCCESS) { return cv_status; }

      auto fn_table = new cvode_user_supplied_fn_table;
      auto cvb_mem = static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(cv_mem, which));
      cvb_mem->python = fn_table;

      cv_status = CVodeSetUserDataB(cv_mem, which, cvb_mem);
      if (cv_status != CV_SUCCESS)
      {
        free(fn_table);
        return cv_status;
      }

      fn_table->fBS = nb::cast(fBS);
      return cv_status;
    },
    nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fBS").none(),
    nb::arg("tB0"), nb::arg("yB0"));

  m.def(
    "CVodeQuadInitBS",
    [](void* cv_mem, int which, std::function<CVQuadRhsStdFnBS> fQBS, N_Vector yQBO)
    {
      auto fn_table = get_cvode_fn_table(cv_mem, which);
      if (fQBS)
      {
        fn_table->fQBS = nb::cast(fQBS);
        return CVodeQuadInitBS(cv_mem, which, cvode_fQBS_wrapper, yQBO);
      }
      else { return CVodeQuadInitBS(cv_mem, which, nullptr, yQBO); }
    },
    nb::arg("cvode_mem"), nb::arg("which"), nb::arg("fQBS").none(),
    nb::arg("yQB0"));

  BIND_CVODEB_CALLBACK(CVodeSetJacFnBS, CVLsJacStdFnBS, lsjacfnBS,
                       cvode_lsjacfnBS_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("jacBS").none());

  BIND_CVODEB_CALLBACK2(CVodeSetPreconditionerBS, CVLsPrecSetupStdFnBS,
                        lsprecsetupfnBS, cvode_lsprecsetupfnBS_wrapper,
                        CVLsPrecSolveStdFnBS, lsprecsolvefnBS,
                        cvode_lsprecsolvefnBS_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("psetBS").none(),
                        nb::arg("psolveBS").none());

  BIND_CVODEB_CALLBACK2(CVodeSetJacTimesBS, CVLsJacTimesSetupStdFnBS,
                        lsjactimessetupfnBS, cvode_lsjactimessetupfnBS_wrapper,
                        CVLsJacTimesVecStdFnBS, lsjactimesvecfnBS,
                        cvode_lsjactimesvecfnBS_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("jsetupBS").none(),
                        nb::arg("jtimesBS").none());

  BIND_CVODEB_CALLBACK(CVodeSetLinSysFnBS, CVLsLinSysStdFnBS, lslinsysfnBS,
                       cvode_lslinsysfnBS_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("linsysBS").none());
}

} // namespace sundials4py

// The destroy functions gets called in our C code by the integrator destructor
extern "C" void cvode_user_supplied_fn_table_destroy(void* ptr)
{
  delete static_cast<cvode_user_supplied_fn_table*>(ptr);
}
