// #ifndef _SUNNONLINSOL_FIXEDPOINT_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNNonlinearSolverContent_FixedPoint =
  nb::class_<_SUNNonlinearSolverContent_FixedPoint>(m, "_SUNNonlinearSolverContent_FixedPoint",
                                                    "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNNonlinSol_FixedPoint",
  [](N_Vector y, int m, SUNContext sunctx)
    -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
  {
    auto SUNNonlinSol_FixedPoint_adapt_return_type_to_shared_ptr =
      [](N_Vector y, int m, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
    {
      auto lambda_result = SUNNonlinSol_FixedPoint(y, m, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNNonlinearSolver>,
                             SUNNonlinearSolverDeleter>(lambda_result);
    };

    return SUNNonlinSol_FixedPoint_adapt_return_type_to_shared_ptr(y, m, sunctx);
  },
  nb::arg("y"), nb::arg("m"), nb::arg("sunctx"), "nb::keep_alive<0, 3>()",
  nb::keep_alive<0, 3>());

m.def(
  "SUNNonlinSol_FixedPointSens",
  [](int count, N_Vector y, int m, SUNContext sunctx)
    -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
  {
    auto SUNNonlinSol_FixedPointSens_adapt_return_type_to_shared_ptr =
      [](int count, N_Vector y, int m, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>>
    {
      auto lambda_result = SUNNonlinSol_FixedPointSens(count, y, m, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNNonlinearSolver>,
                             SUNNonlinearSolverDeleter>(lambda_result);
    };

    return SUNNonlinSol_FixedPointSens_adapt_return_type_to_shared_ptr(count, y,
                                                                       m, sunctx);
  },
  nb::arg("count"), nb::arg("y"), nb::arg("m"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNNonlinSolSetDamping_FixedPoint", SUNNonlinSolSetDamping_FixedPoint,
      nb::arg("NLS"), nb::arg("beta"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
