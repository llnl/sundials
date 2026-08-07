// #ifndef _SUNLINSOL_SPGMR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNLinearSolverContent_SPGMR =
  nb::class_<_SUNLinearSolverContent_SPGMR>(m, "_SUNLinearSolverContent_SPGMR", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNLinSol_SPGMR",
  [](N_Vector y, int pretype, int maxl,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
  {
    auto SUNLinSol_SPGMR_adapt_return_type_to_shared_ptr =
      [](N_Vector y, int pretype, int maxl, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
    {
      auto lambda_result = SUNLinSol_SPGMR(y, pretype, maxl, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNLinearSolver>,
                             SUNLinearSolverDeleter>(lambda_result);
    };

    return SUNLinSol_SPGMR_adapt_return_type_to_shared_ptr(y, pretype, maxl,
                                                           sunctx);
  },
  nb::arg("y"), nb::arg("pretype"), nb::arg("maxl"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNLinSol_SPGMRSetPrecType", SUNLinSol_SPGMRSetPrecType, nb::arg("S"),
      nb::arg("pretype"));

m.def("SUNLinSol_SPGMRSetGSType", SUNLinSol_SPGMRSetGSType, nb::arg("S"),
      nb::arg("gstype"));

m.def("SUNLinSol_SPGMRSetMaxRestarts", SUNLinSol_SPGMRSetMaxRestarts,
      nb::arg("S"), nb::arg("maxrs"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
