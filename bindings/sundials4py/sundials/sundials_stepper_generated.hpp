// #ifndef _SUNDIALS_STEPPER_H
//
// #ifdef __cplusplus
//
// #endif
//

auto pyEnumSUNFullRhsMode = nb::enum_<SUNFullRhsMode>(m, "SUNFullRhsMode",
                                                      nb::is_arithmetic(), "")
                              .value("SUN_FULLRHS_START", SUN_FULLRHS_START, "")
                              .value("SUN_FULLRHS_END", SUN_FULLRHS_END, "")
                              .value("SUN_FULLRHS_OTHER", SUN_FULLRHS_OTHER, "")
                              .export_values();
// #ifndef SWIG
//
// #endif
//

m.def(
  "SUNStepper_Create",
  [](SUNContext sunctx)
    -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
  {
    auto SUNStepper_Create_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNStepper>
    {
      SUNStepper stepper_adapt_modifiable;

      SUNErrCode r = SUNStepper_Create(sunctx, &stepper_adapt_modifiable);
      return std::make_tuple(r, stepper_adapt_modifiable);
    };
    auto SUNStepper_Create_adapt_return_type_to_shared_ptr =
      [&SUNStepper_Create_adapt_modifiable_immutable_to_return](SUNContext sunctx)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
    {
      auto lambda_result =
        SUNStepper_Create_adapt_modifiable_immutable_to_return(sunctx);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNStepper>,
                                             SUNStepperDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SUNStepper_Create_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"),
  "nb::call_policy<sundials4py::returns_references_to<1, 1>>()",
  nb::rv_policy::reference,
  nb::call_policy<sundials4py::returns_references_to<1, 1>>());

m.def(
  "SUNStepper_Evolve",
  [](SUNStepper stepper, sunrealtype tout,
     N_Vector vret) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNStepper_Evolve_adapt_modifiable_immutable_to_return =
      [](SUNStepper stepper, sunrealtype tout,
         N_Vector vret) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      SUNErrCode r = SUNStepper_Evolve(stepper, tout, vret,
                                       &tret_adapt_modifiable);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return SUNStepper_Evolve_adapt_modifiable_immutable_to_return(stepper, tout,
                                                                  vret);
  },
  nb::arg("stepper"), nb::arg("tout"), nb::arg("vret"));

m.def(
  "SUNStepper_OneStep",
  [](SUNStepper stepper, sunrealtype tout,
     N_Vector vret) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNStepper_OneStep_adapt_modifiable_immutable_to_return =
      [](SUNStepper stepper, sunrealtype tout,
         N_Vector vret) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      SUNErrCode r = SUNStepper_OneStep(stepper, tout, vret,
                                        &tret_adapt_modifiable);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return SUNStepper_OneStep_adapt_modifiable_immutable_to_return(stepper,
                                                                   tout, vret);
  },
  nb::arg("stepper"), nb::arg("tout"), nb::arg("vret"));

m.def("SUNStepper_FullRhs", SUNStepper_FullRhs, nb::arg("stepper"),
      nb::arg("t"), nb::arg("v"), nb::arg("f"), nb::arg("mode"));

m.def("SUNStepper_ReInit", SUNStepper_ReInit, nb::arg("stepper"), nb::arg("t0"),
      nb::arg("v0"));

m.def("SUNStepper_Reset", SUNStepper_Reset, nb::arg("stepper"), nb::arg("tR"),
      nb::arg("vR"));

m.def("SUNStepper_ResetCheckpointIndex", SUNStepper_ResetCheckpointIndex,
      nb::arg("stepper"), nb::arg("ckptIdxR"));

m.def("SUNStepper_SetStopTime", SUNStepper_SetStopTime, nb::arg("stepper"),
      nb::arg("tstop"));

m.def("SUNStepper_SetStepDirection", SUNStepper_SetStepDirection,
      nb::arg("stepper"), nb::arg("stepdir"));

m.def(
  "SUNStepper_SetForcing",
  [](SUNStepper stepper, sunrealtype tshift, sunrealtype tscale,
     std::vector<N_Vector> forcing_1d, int nforcing) -> SUNErrCode
  {
    auto SUNStepper_SetForcing_adapt_arr_ptr_to_std_vector =
      [](SUNStepper stepper, sunrealtype tshift, sunrealtype tscale,
         std::vector<N_Vector> forcing_1d, int nforcing) -> SUNErrCode
    {
      N_Vector* forcing_1d_ptr = forcing_1d.empty() ? nullptr : forcing_1d.data();

      auto lambda_result = SUNStepper_SetForcing(stepper, tshift, tscale,
                                                 forcing_1d_ptr, nforcing);
      return lambda_result;
    };

    return SUNStepper_SetForcing_adapt_arr_ptr_to_std_vector(stepper, tshift,
                                                             tscale, forcing_1d,
                                                             nforcing);
  },
  nb::arg("stepper"), nb::arg("tshift"), nb::arg("tscale"),
  nb::arg("forcing_1d"), nb::arg("nforcing"));

m.def("SUNStepper_SetLastFlag", SUNStepper_SetLastFlag, nb::arg("stepper"),
      nb::arg("last_flag"));

m.def(
  "SUNStepper_GetLastFlag",
  [](SUNStepper stepper) -> std::tuple<SUNErrCode, int>
  {
    auto SUNStepper_GetLastFlag_adapt_modifiable_immutable_to_return =
      [](SUNStepper stepper) -> std::tuple<SUNErrCode, int>
    {
      int last_flag_adapt_modifiable;

      SUNErrCode r = SUNStepper_GetLastFlag(stepper, &last_flag_adapt_modifiable);
      return std::make_tuple(r, last_flag_adapt_modifiable);
    };

    return SUNStepper_GetLastFlag_adapt_modifiable_immutable_to_return(stepper);
  },
  nb::arg("stepper"));

m.def(
  "SUNStepper_GetNumSteps",
  [](SUNStepper stepper) -> std::tuple<SUNErrCode, suncountertype>
  {
    auto SUNStepper_GetNumSteps_adapt_modifiable_immutable_to_return =
      [](SUNStepper stepper) -> std::tuple<SUNErrCode, suncountertype>
    {
      suncountertype nst_adapt_modifiable;

      SUNErrCode r = SUNStepper_GetNumSteps(stepper, &nst_adapt_modifiable);
      return std::make_tuple(r, nst_adapt_modifiable);
    };

    return SUNStepper_GetNumSteps_adapt_modifiable_immutable_to_return(stepper);
  },
  nb::arg("stepper"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
