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

#ifndef _SUNDIALS4PY_STEPPER_USERSUPPLIED_HPP
#define _SUNDIALS4PY_STEPPER_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>
#include "sundials4py.hpp"

#include <sundials/sundials_stepper.hpp>

// If helpers are available, include them
#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

struct SUNStepperFunctionTable
{
  nb::object evolve;
  nb::object one_step;
  nb::object full_rhs;
  nb::object reinit;
  nb::object reset;
  nb::object reset_ckpt_idx;
  nb::object set_stop_time;
  nb::object set_step_direction;
  nb::object set_forcing;
  nb::object get_num_steps;
};

template<typename... Args>
inline SUNErrCode sunstepper_evolve_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperEvolveFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::evolve, std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_one_step_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperOneStepFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::one_step, std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_full_rhs_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperFullRhsFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::full_rhs, std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_reinit_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperReInitFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::reinit, std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_reset_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperResetFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::reset, std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_reset_ckpt_idx_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperResetCheckpointIndexFn>,
    SUNStepperFunctionTable, SUNStepper>(&SUNStepperFunctionTable::reset_ckpt_idx,
                                         std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_set_stop_time_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperSetStopTimeFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::set_stop_time,
                std::forward<Args>(args)...);
}

template<typename... Args>
inline SUNErrCode sunstepper_set_step_direction_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNStepperSetStepDirectionFn>, SUNStepperFunctionTable,
    SUNStepper>(&SUNStepperFunctionTable::set_step_direction,
                std::forward<Args>(args)...);
}

using SUNStepperSetForcingStdFn = SUNErrCode(SUNStepper stepper,
                                             sunrealtype tshift,
                                             sunrealtype tscale,
                                             std::vector<N_Vector> forcing,
                                             int nforcing);

inline SUNErrCode sunstepper_set_forcing_wrapper(SUNStepper stepper,
                                                 sunrealtype tshift,
                                                 sunrealtype tscale,
                                                 N_Vector* forcing_1d,
                                                 int nforcing)
{
  auto fn_table = static_cast<SUNStepperFunctionTable*>(stepper->python);
  auto fn =
    nb::cast<std::function<SUNStepperSetForcingStdFn>>(fn_table->set_forcing);

  std::vector<N_Vector> forcing(forcing_1d, forcing_1d + nforcing);

  return fn(stepper, tshift, tscale, forcing, nforcing);
}

using SUNStepperGetNumStepsStdFn =
  std::tuple<SUNErrCode, suncountertype>(SUNStepper);

inline SUNErrCode sunstepper_get_num_steps_wrapper(SUNStepper stepper,
                                                   suncountertype* num_steps)
{
  auto fn_table = static_cast<SUNStepperFunctionTable*>(stepper->python);
  auto fn =
    nb::cast<std::function<SUNStepperGetNumStepsStdFn>>(fn_table->get_num_steps);

  auto result = fn(stepper);

  *num_steps = std::get<1>(result);

  return std::get<0>(result);
}

#endif // _SUNDIALS4PY_STEPPER_USERSUPPLIED_HPP
