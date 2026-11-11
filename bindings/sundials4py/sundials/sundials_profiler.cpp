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
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS SUNProfiler class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"
#include "sundials_comm.hpp"

#include <sundials/sundials_profiler.hpp>

#include "sundials/sundials_types.h"
#include "sundials_profiler_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunprofiler(nb::module_& m)
{
#include "sundials_profiler_generated.hpp"

  m.def(
    "SUNProfiler_Create",
    [](nb::handle comm, const char* title)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNProfiler>>>
    {
      SUNComm c = to_suncomm(comm);
      SUNProfiler p_out;
      SUNErrCode r = SUNProfiler_Create(c, title, &p_out);
      return std::make_tuple(r,
                             our_make_shared<std::remove_pointer_t<SUNProfiler>,
                                             SUNProfilerDeleter>(p_out));
    },
    nb::arg("comm").none(), nb::arg("title"), nb::rv_policy::reference,
    "Create a SUNProfiler. With a non-null communicator the call is "
    "collective over comm (the profiler duplicates it).");

  nb::class_<SUNProfiler_>(m, "SUNProfiler_");
}

} // namespace sundials4py
