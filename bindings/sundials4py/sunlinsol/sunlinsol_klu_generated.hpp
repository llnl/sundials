// #ifndef _SUNLINSOL_KLU_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("SUNKLU_ORDERING_DEFAULT") = SUNKLU_ORDERING_DEFAULT;

m.attr("SUNKLU_REINIT_FULL") = SUNKLU_REINIT_FULL;

m.attr("SUNKLU_REINIT_PARTIAL") = SUNKLU_REINIT_PARTIAL;

auto pyClass_SUNLinearSolverContent_KLU =
  nb::class_<_SUNLinearSolverContent_KLU>(m, "_SUNLinearSolverContent_KLU", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNLinSol_KLU",
  [](N_Vector y, SUNMatrix A,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
  {
    auto SUNLinSol_KLU_adapt_return_type_to_shared_ptr =
      [](N_Vector y, SUNMatrix A, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
    {
      auto lambda_result = SUNLinSol_KLU(y, A, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNLinearSolver>,
                             SUNLinearSolverDeleter>(lambda_result);
    };

    return SUNLinSol_KLU_adapt_return_type_to_shared_ptr(y, A, sunctx);
  },
  nb::arg("y"), nb::arg("A"), nb::arg("sunctx"), "nb::keep_alive<0, 3>()",
  nb::keep_alive<0, 3>());

m.def("SUNLinSol_KLUReInit", SUNLinSol_KLUReInit, nb::arg("S"), nb::arg("A"),
      nb::arg("nnz"), nb::arg("reinit_type"));

m.def("SUNLinSol_KLUSetOrdering", SUNLinSol_KLUSetOrdering, nb::arg("S"),
      nb::arg("ordering_choice"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
