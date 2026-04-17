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
 * SUNDIALS SUNLogger class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_logger.hpp>
#include <sundials/sundials_types.h>

#include "sundials_logger_impl.h"
#include "sundials_logger_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunlogger(nb::module_& m)
{
#include "sundials_logger_generated.hpp"
  nb::class_<SUNLogger_>(m, "SUNLogger_");

  m.def("SUNLogger_SetErrorFile", SUNLogger_SetErrorFile, nb::arg("logger"),
        nb::arg("error_fp").none());

  m.def("SUNLogger_SetWarningFile", SUNLogger_SetWarningFile, nb::arg("logger"),
        nb::arg("warning_fp").none());

  m.def("SUNLogger_SetDebugFile", SUNLogger_SetDebugFile, nb::arg("logger"),
        nb::arg("debug_fp").none());

  m.def("SUNLogger_SetInfoFile", SUNLogger_SetInfoFile, nb::arg("logger"),
        nb::arg("info_fp").none());

  m.def(
    "SUNLogger_SetQueueAndFlushMsgFns",
    [](SUNLogger logger,
       std::function<std::remove_pointer_t<SUNLoggerQueueMsgFn>> queue_fn,
       std::function<std::remove_pointer_t<SUNLoggerFlushMsgFn>> flush_fn) -> SUNErrCode
    {
      if (!logger->python) { logger->python = new SUNLoggerFunctionTable; }
      auto fntable = static_cast<SUNLoggerFunctionTable*>(logger->python);
      fntable->queue_msg = nb::cast(queue_fn);
      fntable->flush_msg = nb::cast(flush_fn);
      if (queue_fn && flush_fn)
      {
        return SUNLogger_SetQueueAndFlushMsgFns(logger,
                                                sunlogger_queue_msg_wrapper,
                                                sunlogger_flush_msg_wrapper,
                                                nullptr);
      }
      else if (queue_fn && !flush_fn)
      {
        return SUNLogger_SetQueueAndFlushMsgFns(logger,
                                                sunlogger_queue_msg_wrapper,
                                                nullptr, nullptr);
      }
      else
      {
        return SUNLogger_SetQueueAndFlushMsgFns(logger, nullptr, nullptr, nullptr);
      }
    },
    nb::arg("logger"), nb::arg("queue_fn").none(), nb::arg("flush_fn").none());
}

} // namespace sundials4py

extern "C" void SUNLoggerFunctionTable_Destroy(void* ptr)
{
  delete static_cast<SUNLoggerFunctionTable*>(ptr);
}

