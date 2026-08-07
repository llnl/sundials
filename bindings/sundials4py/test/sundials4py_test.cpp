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
 * This file defines the sundials4py test module, which contains
 * utility functions for testing purposes.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_context.hpp>
#include <sundials/sundials_errors.h>

#include <sundials/priv/sundials_errors_impl.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_test(nb::module_& m)
{
  m.def(
    "SUNContext_TestErrHandler",
    [](SUNContext sunctx)
    {
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__,
                          "create an error to test the error handlers",
                          SUN_ERR_ARG_CORRUPT, sunctx);
    },
    "This function is for testing purposes and should not be called.");
}

} // namespace sundials4py
