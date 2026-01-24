// #ifndef _SUNLINSOL_SPTFQMR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNLinearSolverContent_SPTFQMR =
  nb::class_<_SUNLinearSolverContent_SPTFQMR>(m,
                                              "_SUNLinearSolverContent_SPTFQMR",
                                              "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNLinSol_SPTFQMR",
  [](N_Vector y, int pretype, int maxl,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
  {
    auto SUNLinSol_SPTFQMR_adapt_return_type_to_shared_ptr =
      [](N_Vector y, int pretype, int maxl, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>>
    {
      auto lambda_result = SUNLinSol_SPTFQMR(y, pretype, maxl, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNLinearSolver>,
                             SUNLinearSolverDeleter>(lambda_result);
    };

    return SUNLinSol_SPTFQMR_adapt_return_type_to_shared_ptr(y, pretype, maxl,
                                                             sunctx);
  },
  nb::arg("y"), nb::arg("pretype"), nb::arg("maxl"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNLinSol_SPTFQMRSetPrecType", SUNLinSol_SPTFQMRSetPrecType,
      nb::arg("S"), nb::arg("pretype"));

m.def("SUNLinSol_SPTFQMRSetMaxl", SUNLinSol_SPTFQMRSetMaxl, nb::arg("S"),
      nb::arg("maxl"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
