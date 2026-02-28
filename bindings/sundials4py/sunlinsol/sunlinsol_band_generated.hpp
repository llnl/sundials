// #ifndef _SUNLINSOL_BAND_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNLinearSolverContent_Band =
  nb::class_<_SUNLinearSolverContent_Band>(m, "_SUNLinearSolverContent_Band", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNLinSol_Band",
  [](N_Vector y, SUNMatrix A,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
  {
    auto SUNLinSol_Band_adapt_return_type_to_shared_ptr =
      [](N_Vector y, SUNMatrix A, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
    {
      auto lambda_result = SUNLinSol_Band(y, A, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNLinearSolver>,
                             SUNLinearSolverDeleter>(lambda_result);
    };

    return SUNLinSol_Band_adapt_return_type_to_shared_ptr(y, A, sunctx);
  },
  nb::arg("y"), nb::arg("A"), nb::arg("sunctx"), "nb::keep_alive<0, 3>()",
  nb::keep_alive<0, 3>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
