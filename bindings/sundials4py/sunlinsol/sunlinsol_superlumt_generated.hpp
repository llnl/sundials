// #ifndef _SUNLINSOL_SLUMT_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNLinearSolverContent_SuperLUMT =
  nb::class_<_SUNLinearSolverContent_SuperLUMT>(m, "_SUNLinearSolverContent_SuperLUMT",
                                                "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNLinSol_SuperLUMT",
  [](N_Vector y, SUNMatrix A, int num_threads,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
  {
    auto SUNLinSol_SuperLUMT_adapt_return_type_to_shared_ptr =
      [](N_Vector y, SUNMatrix A, int num_threads, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
    {
      auto lambda_result = SUNLinSol_SuperLUMT(y, A, num_threads, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNLinearSolver>,
                             SUNLinearSolverDeleter>(lambda_result);
    };

    return SUNLinSol_SuperLUMT_adapt_return_type_to_shared_ptr(y, A, num_threads,
                                                               sunctx);
  },
  nb::arg("y"), nb::arg("A"), nb::arg("num_threads"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNLinSol_SuperLUMTSetOrdering", SUNLinSol_SuperLUMTSetOrdering,
      nb::arg("S"), nb::arg("ordering_choice"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
