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
 * SUNDIALS SUNStepper class. It contains hand-written code 
 * for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_stepper.hpp>
#include <sundials/sundials_types.h>

#include "sundials_stepper_impl.h"
#include "sundials_stepper_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunstepper(nb::module_& m)
{
#include "sundials_stepper_generated.hpp"

  nb::class_<SUNStepper_>(m, "SUNStepper_");

  m.def(
    "SUNStepper_SetEvolveFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperEvolveFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable    = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->evolve = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetEvolveFn(stepper, sunstepper_evolve_wrapper);
      }
      else { return SUNStepper_SetEvolveFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetOneStepFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperOneStepFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->one_step = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetOneStepFn(stepper, sunstepper_one_step_wrapper);
      }
      else { return SUNStepper_SetOneStepFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetFullRhsFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperFullRhsFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->full_rhs = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetFullRhsFn(stepper, sunstepper_full_rhs_wrapper);
      }
      else { return SUNStepper_SetFullRhsFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetReInitFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperReInitFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable    = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->reinit = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetReInitFn(stepper, sunstepper_reinit_wrapper);
      }
      else { return SUNStepper_SetReInitFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetResetFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperResetFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable   = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->reset = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetResetFn(stepper, sunstepper_reset_wrapper);
      }
      else { return SUNStepper_SetResetFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetResetCheckpointIndexFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperResetCheckpointIndexFn>> fn)
      -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->reset_ckpt_idx = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetResetCheckpointIndexFn(stepper,
                                                    sunstepper_reset_ckpt_idx_wrapper);
      }
      else { return SUNStepper_SetResetCheckpointIndexFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetStopTimeFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperSetStopTimeFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->set_stop_time = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetStopTimeFn(stepper,
                                        sunstepper_set_stop_time_wrapper);
      }
      else { return SUNStepper_SetStopTimeFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetStepDirectionFn",
    [](SUNStepper stepper,
       std::function<std::remove_pointer_t<SUNStepperSetStepDirectionFn>> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->set_step_direction = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetStepDirectionFn(stepper,
                                             sunstepper_set_step_direction_wrapper);
      }
      else { return SUNStepper_SetStepDirectionFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetForcingFn",
    [](SUNStepper stepper, std::function<SUNStepperSetForcingStdFn> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->set_forcing = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetForcingFn(stepper, sunstepper_set_forcing_wrapper);
      }
      else { return SUNStepper_SetForcingFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());

  m.def(
    "SUNStepper_SetGetNumStepsFn",
    [](SUNStepper stepper, std::function<SUNStepperGetNumStepsStdFn> fn) -> SUNErrCode
    {
      if (!stepper->python) { stepper->python = new SUNStepperFunctionTable; }
      auto fntable = static_cast<SUNStepperFunctionTable*>(stepper->python);
      fntable->get_num_steps = nb::cast(fn);
      if (fn)
      {
        return SUNStepper_SetGetNumStepsFn(stepper,
                                           sunstepper_get_num_steps_wrapper);
      }
      else { return SUNStepper_SetGetNumStepsFn(stepper, nullptr); }
    },
    nb::arg("stepper"), nb::arg("fn").none());
}

} // namespace sundials4py

extern "C" void SUNStepperFunctionTable_Destroy(void* ptr)
{
  delete static_cast<SUNStepperFunctionTable*>(ptr);
}
