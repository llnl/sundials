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

auto pyClass_SUNNonlinearSolverContent_Auto =
  nb::class_<_SUNNonlinearSolverContent_Auto>(m,
                                              "_SUNNonlinearSolverContent_Auto",
                                              "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNNonlinSol_Auto",
  [](N_Vector y, int m, SUNNonlinSolAutoType type, SUNContext sunctx)
    -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
  {
    auto SUNNonlinSol_Auto_adapt_return_type_to_shared_ptr =
      [](N_Vector y, int m, SUNNonlinSolAutoType type, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
    {
      auto lambda_result = SUNNonlinSol_Auto(y, m, type, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNNonlinearSolver>,
                             SUNNonlinearSolverDeleter>(lambda_result);
    };

    return SUNNonlinSol_Auto_adapt_return_type_to_shared_ptr(y, m, type, sunctx);
  },
  nb::arg("y"), nb::arg("m"), nb::arg("type"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
