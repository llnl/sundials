// #ifndef _SUNNONLINSOL_NEWTON_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNNonlinearSolverContent_Newton =
  nb::class_<_SUNNonlinearSolverContent_Newton>(m, "_SUNNonlinearSolverContent_Newton",
                                                "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNNonlinSol_Newton",
  [](N_Vector y, SUNContext sunctx)
    -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
  {
    auto SUNNonlinSol_Newton_adapt_return_type_to_shared_ptr =
      [](N_Vector y, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
    {
      auto lambda_result = SUNNonlinSol_Newton(y, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNNonlinearSolver>,
                             SUNNonlinearSolverDeleter>(lambda_result);
    };

    return SUNNonlinSol_Newton_adapt_return_type_to_shared_ptr(y, sunctx);
  },
  nb::arg("y"), nb::arg("sunctx"), "nb::keep_alive<0, 2>()",
  nb::keep_alive<0, 2>());

m.def(
  "SUNNonlinSol_NewtonSens",
  [](int count, N_Vector y, SUNContext sunctx)
    -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
  {
    auto SUNNonlinSol_NewtonSens_adapt_return_type_to_shared_ptr =
      [](int count, N_Vector y, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
    {
      auto lambda_result = SUNNonlinSol_NewtonSens(count, y, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNNonlinearSolver>,
                             SUNNonlinearSolverDeleter>(lambda_result);
    };

    return SUNNonlinSol_NewtonSens_adapt_return_type_to_shared_ptr(count, y,
                                                                   sunctx);
  },
  nb::arg("count"), nb::arg("y"), nb::arg("sunctx"), "nb::keep_alive<0, 3>()",
  nb::keep_alive<0, 3>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
