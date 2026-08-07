// #ifndef _SUNLINSOL_SPFGMR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNLinearSolverContent_SPFGMR =
  nb::class_<_SUNLinearSolverContent_SPFGMR>(m, "_SUNLinearSolverContent_SPFGMR",
                                             "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNLinSol_SPFGMR",
  [](N_Vector y, int pretype, int maxl,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
  {
    auto SUNLinSol_SPFGMR_adapt_return_type_to_shared_ptr =
      [](N_Vector y, int pretype, int maxl, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
    {
      auto lambda_result = SUNLinSol_SPFGMR(y, pretype, maxl, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNLinearSolver>,
                             SUNLinearSolverDeleter>(lambda_result);
    };

    return SUNLinSol_SPFGMR_adapt_return_type_to_shared_ptr(y, pretype, maxl,
                                                            sunctx);
  },
  nb::arg("y"), nb::arg("pretype"), nb::arg("maxl"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNLinSol_SPFGMRSetPrecType", SUNLinSol_SPFGMRSetPrecType, nb::arg("S"),
      nb::arg("pretype"));

m.def("SUNLinSol_SPFGMRSetGSType", SUNLinSol_SPFGMRSetGSType, nb::arg("S"),
      nb::arg("gstype"));

m.def("SUNLinSol_SPFGMRSetMaxRestarts", SUNLinSol_SPFGMRSetMaxRestarts,
      nb::arg("S"), nb::arg("maxrs"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
