/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_LOGGER_USERSUPPLIED_HPP
#define _SUNDIALS4PY_LOGGER_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>
#include "sundials4py.hpp"

#include <sundials/sundials_logger.h>

// If helpers are available, include them
#include "sundials4py_helpers.hpp"

namespace nb = nanobind;

struct SUNLoggerFunctionTable
{
  nb::object queue_msg;
  nb::object flush_msg;
};

template<typename... Args>
inline SUNErrCode sunlogger_queue_msg_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNLoggerQueueMsgFn>, SUNLoggerFunctionTable,
    SUNLogger>(&SUNLoggerFunctionTable::queue_msg, std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunlogger_flush_msg_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNLoggerFlushMsgFn>, SUNLoggerFunctionTable,
    SUNLogger>(&SUNLoggerFunctionTable::flush_msg, std::forward<Args>(args)...);
}

#endif // _SUNDIALS4PY_LOGGER_USERSUPPLIED_HPP
