/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This file defines the sundials4py.core module. 
 * -----------------------------------------------------------------*/

#include "sundials/sundials_types.h"
#include "sundials4py.hpp"

#include <stdexcept>
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_futils.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector(nb::module_& m);
void bind_sunadaptcontroller(nb::module_& m);
void bind_sunadjointcheckpointscheme(nb::module_& m);
void bind_sunadjointstepper(nb::module_& m);
void bind_suncontext(nb::module_& m);
void bind_sundomeigestimator(nb::module_& m);
void bind_sunlinearsolver(nb::module_& m);
void bind_sunlogger(nb::module_& m);
void bind_sunmatrix(nb::module_& m);
void bind_sunmemory(nb::module_& m);
void bind_sunnonlinearsolver(nb::module_& m);
void bind_sunprofiler(nb::module_& m);
void bind_sunstepper(nb::module_& m);

void bind_core(nb::module_& m)
{
#include "sundials_errors.hpp"
#include "sundials_types_generated.hpp"

  m.attr("SUNDIALS_VERSION")     = SUNDIALS_VERSION;
  m.attr("SUNDIALS_GIT_VERSION") = SUNDIALS_GIT_VERSION;

  // handle opening and closing C files
  nb::class_<FILE>(m, "FILE");
  m.def("SUNFileOpen",
        [](const char* filename, const char* modes)
        {
          FILE* tmp = nullptr;
          std::shared_ptr<FILE> fp;
          SUNErrCode status = SUNFileOpen(filename, modes, &tmp);
          if (status) { fp = nullptr; }
          else { fp = std::shared_ptr<FILE>(tmp, std::fclose); }
          return std::make_tuple(status, fp);
        });

  bind_nvector(m);
  bind_sunadaptcontroller(m);
  bind_sunadjointcheckpointscheme(m);
  bind_sunadjointstepper(m);
  bind_suncontext(m);
  bind_sundomeigestimator(m);
  bind_sunlinearsolver(m);
  bind_sunlogger(m);
  bind_sunmatrix(m);
  bind_sunmemory(m);
  bind_sunnonlinearsolver(m);
  bind_sunprofiler(m);
  bind_sunstepper(m);

  //
  // Expose sunrealtye and sunindextype as the corresponding numpy types
  //

  nb::object np = nb::module_::import_("numpy");
#if defined(SUNDIALS_SINGLE_PRECISION)
  m.attr("sunrealtype") = np.attr("float32");
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  m.attr("sunrealtype") = np.attr("float64");
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  m.attr("sunrealtype") = np.attr("longdouble");
#else
#error Unknown sunrealtype, email sundials-users@llnl.gov
#endif

#if defined(SUNDIALS_INT64_T)
  m.attr("sunindextype") = np.attr("int64");
#elif defined(SUNDIALS_INT32_T)
  m.attr("sunindextype") = np.attr("int32");
#else
#error Unknown sunindextype, email sundials-users@llnl.gov
#endif

  if (sizeof(suncountertype) == 8)
  {
    m.attr("suncountertype") = np.attr("int64");
  }
  else if (sizeof(suncountertype) == 4)
  {
    m.attr("suncountertype") = np.attr("int32");
  }
  else
  {
    throw sundials4py::illegal_value(
      "unknown suncountertype,  email sundials-users@llnl.gov\n");
  }
}

} // namespace sundials4py
