// #ifndef SUNDIALS_NONLINSOLAUTO_H_
//
// #ifdef __cplusplus
//
// #endif
//

auto pyEnumSUNNonlinSolAutoType =
  nb::enum_<SUNNonlinSolAutoType>(m, "SUNNonlinSolAutoType",
                                  nb::is_arithmetic(), "")
    .value("SUNNONLINSOL_AUTO_FIXEDPOINT", SUNNONLINSOL_AUTO_FIXEDPOINT, "")
    .value("SUNNONLINSOL_AUTO_NEWTON", SUNNONLINSOL_AUTO_NEWTON, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClassSUNNonlinearSolverContent_Auto_ =
  nb::class_<SUNNonlinearSolverContent_Auto_>(m,
                                              "SUNNonlinearSolverContent_Auto_",
                                              "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNNonlinSol_Auto",
  [](N_Vector y, int m, SUNNonlinSolAutoType initial_solver_type,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
  {
    auto SUNNonlinSol_Auto_adapt_return_type_to_shared_ptr =
      [](N_Vector y, int m, SUNNonlinSolAutoType initial_solver_type,
         SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
    {
      auto lambda_result = SUNNonlinSol_Auto(y, m, initial_solver_type, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNNonlinearSolver>,
                             SUNNonlinearSolverDeleter>(lambda_result);
    };

    return SUNNonlinSol_Auto_adapt_return_type_to_shared_ptr(y, m,
                                                             initial_solver_type,
                                                             sunctx);
  },
  nb::arg("y"), nb::arg("m"), nb::arg("initial_solver_type"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNNonlinSolSetSwitchingParameters_Auto",
      SUNNonlinSolSetSwitchingParameters_Auto, nb::arg("NLS"),
      nb::arg("newt_to_fp_threshold"), nb::arg("newt_to_fp_delay"),
      nb::arg("fp_to_newt_threshold"), nb::arg("fp_to_newt_delay"));

m.def(
  "SUNNonlinSolGetNumItersByType_Auto",
  [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long, long>
  {
    auto SUNNonlinSolGetNumItersByType_Auto_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long, long>
    {
      long fp_iters_adapt_modifiable;
      long newt_iters_adapt_modifiable;

      SUNErrCode r =
        SUNNonlinSolGetNumItersByType_Auto(NLS, &fp_iters_adapt_modifiable,
                                           &newt_iters_adapt_modifiable);
      return std::make_tuple(r, fp_iters_adapt_modifiable,
                             newt_iters_adapt_modifiable);
    };

    return SUNNonlinSolGetNumItersByType_Auto_adapt_modifiable_immutable_to_return(
      NLS);
  },
  nb::arg("NLS"));

m.def(
  "SUNNonlinSolGetNumConvFailsByType_Auto",
  [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long, long>
  {
    auto SUNNonlinSolGetNumConvFailsByType_Auto_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long, long>
    {
      long fp_nconvfails_adapt_modifiable;
      long newt_nconvfails_adapt_modifiable;

      SUNErrCode r =
        SUNNonlinSolGetNumConvFailsByType_Auto(NLS,
                                               &fp_nconvfails_adapt_modifiable,
                                               &newt_nconvfails_adapt_modifiable);
      return std::make_tuple(r, fp_nconvfails_adapt_modifiable,
                             newt_nconvfails_adapt_modifiable);
    };

    return SUNNonlinSolGetNumConvFailsByType_Auto_adapt_modifiable_immutable_to_return(
      NLS);
  },
  nb::arg("NLS"));

m.def(
  "SUNNonlinSolGetUpdateNorm_Auto",
  [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNNonlinSolGetUpdateNorm_Auto_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype delnrm_adapt_modifiable;

      SUNErrCode r = SUNNonlinSolGetUpdateNorm_Auto(NLS,
                                                    &delnrm_adapt_modifiable);
      return std::make_tuple(r, delnrm_adapt_modifiable);
    };

    return SUNNonlinSolGetUpdateNorm_Auto_adapt_modifiable_immutable_to_return(
      NLS);
  },
  nb::arg("NLS"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
