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

#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointcheckpointscheme.hpp>
#include <sundials/sundials_core.hpp>

#include "sundials/sundials_adjointcheckpointscheme.h"
#include "sundials_adjointcheckpointscheme_impl.h"

#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunadjointcheckpointscheme_fixed(nb::module_& m)
{
#include "sunadjointcheckpointscheme_fixed_generated.hpp"

  m.def(
    "SUNAdjointCheckpointScheme_Create_Fixed",
    [](SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval,
       suncountertype estimate, sunbooleantype keep, SUNContext sunctx)
      -> std::tuple<SUNErrCode,
                    std::shared_ptr<std::remove_pointer_t<SUNAdjointCheckpointScheme>>>
    {
      SUNAdjointCheckpointScheme check_scheme;
      SUNErrCode err =
        SUNAdjointCheckpointScheme_Create_Fixed(io_mode, mem_helper, interval,
                                                estimate, keep, sunctx,
                                                &check_scheme);
      return std::make_tuple(err,
                             our_make_shared<
                               std::remove_pointer_t<SUNAdjointCheckpointScheme>,
                               SUNAdjointCheckpointSchemeDeleter>(check_scheme));
    },
    nb::arg("io_mode"), nb::arg("mem_helper"), nb::arg("interval"),
    nb::arg("estimate"), nb::arg("keep"), nb::arg("sunctx"),
    nb::call_policy<sundials4py::returns_references_to<
      2, 1>>() /* keep SUNMemoryHelper alive as long as SUNAdjointCheckpointScheme is alive */,
    nb::call_policy<sundials4py::returns_references_to<
      6, 1>>() /* keep SUNContext alive as long as SUNAdjointCheckpointScheme is alive */);
}

} // namespace sundials4py
