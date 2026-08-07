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
 * SUNDIALS SUNAdjointCheckpointScheme class. It contains hand-written
 * code for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_adjointcheckpointscheme.hpp>
#include <sundials/sundials_nvector.hpp>

#include "sundials_adjointcheckpointscheme_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunadjointcheckpointscheme(nb::module_& m)
{
#include "sundials_adjointcheckpointscheme_generated.hpp"

  nb::class_<SUNAdjointCheckpointScheme_>(m, "SUNAdjointCheckpointScheme_");

  m.def(
    "SUNAdjointCheckpointScheme_LoadVector",
    [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
       suncountertype stage_num, sunbooleantype peek, N_Vector tmpl)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<N_Vector>>, sunrealtype>
    {
      N_Vector nvec_out = N_VClone(tmpl);
      sunrealtype tout;
      SUNErrCode err =
        SUNAdjointCheckpointScheme_LoadVector(check_scheme, step_num, stage_num,
                                              peek, &nvec_out, &tout);

      return std::make_tuple(err,
                             our_make_shared<std::remove_pointer_t<N_Vector>,
                                             N_VectorDeleter>(nvec_out),
                             tout);
    },
    nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"),
    nb::arg("peek"), nb::arg("tmpl"));
}

} // namespace sundials4py
