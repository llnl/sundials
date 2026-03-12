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
 * SUNDIALS SUNContext class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials/sundials_errors.h"
#include "sundials4py.hpp"

#include <sundials/sundials_context.hpp>
#include <sundials/sundials_logger.hpp>
#include <sundials/sundials_profiler.hpp>

#include "sundials/sundials_types.h"
#include "sundials_logger_impl.h"
#include "sundials_profiler_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;
using sundials::SUNContextDeleter;

#include "sundials_context_usersupplied.hpp"

namespace sundials4py {

void bind_suncontext(nb::module_& m)
{
#include "sundials_context_generated.hpp"

  nb::class_<SUNContext_>(m, "SUNContext_");

  m.def("SUNContext_PushErrHandler",
        [](SUNContext sunctx,
           std::function<std::remove_pointer_t<SUNErrHandlerFn>> err_fn)
        {
          if (!err_fn) { throw sundials4py::illegal_value("err_fn was None"); }

          if (!sunctx->python)
          {
            sunctx->python = new SUNContextFunctionTable;

            // Only push the wrapper the first time this is called
            SUNErrCode status =
              SUNContext_PushErrHandler(sunctx, suncontext_errhandler_wrapper,
                                        sunctx->python);
            if (status)
            {
              throw sundials4py::error_returned(
                "SUNContext_PushErrHandler returned an error");
            }
          }

          auto fn_table = static_cast<SUNContextFunctionTable*>(sunctx->python);

          fn_table->err_handlers.push_back(nb::cast(err_fn));

          return SUN_SUCCESS;
        });

  m.def("SUNContext_PopErrHandler",
        [](SUNContext sunctx) -> SUNErrCode
        {
          if (!sunctx->python) { return SUN_SUCCESS; }

          auto fn_table = static_cast<SUNContextFunctionTable*>(sunctx->python);

          if (fn_table->err_handlers.size() > 0)
          {
            // pop the python functions off the interface layer stack
            fn_table->err_handlers.pop_back();
          }

          if (fn_table->err_handlers.size() == 0)
          {
            // now we can pop the suncontext_errhandler_wrapper off the C side stack
            return SUNContext_PopErrHandler(sunctx);
          }

          return SUN_SUCCESS;
        });
}

} // namespace sundials4py

extern "C" void SUNContextFunctionTable_Destroy(void* ptr)
{
  delete static_cast<SUNContextFunctionTable*>(ptr);
}
