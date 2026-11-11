/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This file defines the sundials4py Python module and includes all
 * of the submodule pieces.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"
#include "sundials/sundials_config.h"

namespace nb = nanobind;

namespace sundials4py {

//
// Forward declarations of all of the binding functions
//

void bind_core(nb::module_& m);

void bind_test(nb::module_& m);

void bind_arkode(nb::module_& m);
void bind_cvodes(nb::module_& m);
void bind_idas(nb::module_& m);
void bind_kinsol(nb::module_& m);

void bind_nvector_serial(nb::module_& m);
void bind_nvector_manyvector(nb::module_& m);
#ifdef SUNDIALS_NVECTOR_MPIMANYVECTOR
void bind_nvector_mpimanyvector(nb::module_& m);
#endif
#ifdef SUNDIALS_NVECTOR_MPIPLUSX
void bind_nvector_mpiplusx(nb::module_& m);
#endif
#ifdef SUNDIALS4PY_NVECTOR_CUDA
void bind_nvector_cuda(nb::module_& m);
#endif
#ifdef SUNDIALS4PY_NVECTOR_HIP
void bind_nvector_hip(nb::module_& m);
#endif

void bind_sumemoryhelper_sys(nb::module_& m);
#ifdef SUNDIALS4PY_NVECTOR_CUDA
void bind_sunmemoryhelper_cuda(nb::module_& m);
#endif
#ifdef SUNDIALS4PY_NVECTOR_HIP
void bind_sunmemoryhelper_hip(nb::module_& m);
#endif

void bind_sunadaptcontroller_imexgus(nb::module_& m);
void bind_sunadaptcontroller_mrihtol(nb::module_& m);
void bind_sunadaptcontroller_soderlind(nb::module_& m);

void bind_sunadjointcheckpointscheme_fixed(nb::module_& m);

void bind_sundomeigest_power(nb::module_& m);

void bind_sunlinsol_band(nb::module_& m);
void bind_sunlinsol_dense(nb::module_& m);
#if defined(SUNDIALS_KLU_ENABLED)
void bind_sunlinsol_klu(nb::module_& m);
#endif
void bind_sunlinsol_pcg(nb::module_& m);
void bind_sunlinsol_spbcgs(nb::module_& m);
void bind_sunlinsol_spfgmr(nb::module_& m);
void bind_sunlinsol_spgmr(nb::module_& m);
void bind_sunlinsol_sptfqmr(nb::module_& m);

void bind_sunmatrix_band(nb::module_& m);
void bind_sunmatrix_dense(nb::module_& m);
void bind_sunmatrix_sparse(nb::module_& m);

void bind_sunnonlinsol_fixedpoint(nb::module_& m);
void bind_sunnonlinsol_newton(nb::module_& m);
void bind_sunnonlinsol_auto(nb::module_& m);

} // namespace sundials4py

//
// Define main module, sundials4py, and all of its submodules
//

#ifdef SUNDIALS4PY_ROCM_DISTRIBUTION
NB_MODULE(_native, m)
#else
NB_MODULE(sundials4py, m)
#endif
{
#ifndef SUNDIALS4PY_ROCM_DISTRIBUTION
  // The optional ROCm wheel owns a separate extension. Load it before any
  // CPU bindings are registered so its nanobind types are the only active
  // SUNDIALS types in this process.
  nb::object rocm_spec =
    nb::module_::import_("importlib.util").attr("find_spec")(
      "sundials4py_rocm");
  if (!rocm_spec.is_none())
  {
    nb::module_ rocm = nb::module_::import_("sundials4py_rocm");
    nb::dict modules = nb::module_::import_("sys").attr("modules");
    m.attr("__version__") = rocm.attr("__version__");
    for (const char* name : {"core", "test", "arkode", "cvodes", "idas",
                             "kinsol"})
    {
      nb::object submodule = rocm.attr(name);
      m.attr(name)         = submodule;
      modules[nb::str((std::string("sundials4py.") + name).c_str())] =
        submodule;
    }
    return;
  }
#endif
#ifdef NDEBUG
  // The nanobind leak warnings can be quite noisy due to leaks within Python itself, so we disable them for Release builds.
  nb::set_leak_warnings(false);
#endif

  // Setting the version in __version__ is a quasi-standard practice for Python modules
  m.attr("__version__") = SUNDIALS_VERSION;

  nb::module_ core_m = m.def_submodule("core", "A submodule of 'sundials4py'");
  sundials4py::bind_core(core_m);

  //
  // Create test submodule
  //

  nb::module_ test_m =
    m.def_submodule("test", "A submodule of 'sundials4py' for testing");
  sundials4py::bind_test(test_m);

  //
  // Create submodules for each package
  //

  nb::module_ arkode_m = m.def_submodule("arkode",
                                         "A submodule of 'sundials4py'");
  sundials4py::bind_arkode(arkode_m);

  nb::module_ cvodes_m = m.def_submodule("cvodes",
                                         "A submodule of 'sundials4py'");
  sundials4py::bind_cvodes(cvodes_m);

  nb::module_ idas_m = m.def_submodule("idas", "A submodule of 'sundials4py'");
  sundials4py::bind_idas(idas_m);

  nb::module_ kinsol_m = m.def_submodule("kinsol",
                                         "A submodule of 'sundials4py'");
  sundials4py::bind_kinsol(kinsol_m);

  //
  // Bind all implementation modules directly to core_m
  //

  sundials4py::bind_nvector_serial(core_m);
  sundials4py::bind_nvector_manyvector(core_m);
#ifdef SUNDIALS_NVECTOR_MPIMANYVECTOR
  sundials4py::bind_nvector_mpimanyvector(core_m);
#endif
#ifdef SUNDIALS_NVECTOR_MPIPLUSX
  sundials4py::bind_nvector_mpiplusx(core_m);
#endif
#ifdef SUNDIALS4PY_NVECTOR_CUDA
  sundials4py::bind_nvector_cuda(core_m);
#endif
#ifdef SUNDIALS4PY_NVECTOR_HIP
  sundials4py::bind_nvector_hip(core_m);
#endif

  sundials4py::bind_sunadaptcontroller_imexgus(core_m);
  sundials4py::bind_sunadaptcontroller_mrihtol(core_m);
  sundials4py::bind_sunadaptcontroller_soderlind(core_m);

  sundials4py::bind_sunadjointcheckpointscheme_fixed(core_m);

  sundials4py::bind_sundomeigest_power(core_m);

  sundials4py::bind_sunlinsol_band(core_m);
  sundials4py::bind_sunlinsol_dense(core_m);
#if defined(SUNDIALS_KLU_ENABLED)
  sundials4py::bind_sunlinsol_klu(core_m);
#endif
  sundials4py::bind_sunlinsol_pcg(core_m);
  sundials4py::bind_sunlinsol_spbcgs(core_m);
  sundials4py::bind_sunlinsol_spfgmr(core_m);
  sundials4py::bind_sunlinsol_spgmr(core_m);
  sundials4py::bind_sunlinsol_sptfqmr(core_m);

  sundials4py::bind_sunmatrix_band(core_m);
  sundials4py::bind_sunmatrix_dense(core_m);
  sundials4py::bind_sunmatrix_sparse(core_m);

  sundials4py::bind_sumemoryhelper_sys(core_m);
#ifdef SUNDIALS4PY_NVECTOR_CUDA
  sundials4py::bind_sunmemoryhelper_cuda(core_m);
#endif
#ifdef SUNDIALS4PY_NVECTOR_HIP
  sundials4py::bind_sunmemoryhelper_hip(core_m);
#endif

  sundials4py::bind_sunnonlinsol_fixedpoint(core_m);
  sundials4py::bind_sunnonlinsol_newton(core_m);
  sundials4py::bind_sunnonlinsol_auto(core_m);
}
