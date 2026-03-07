// #ifndef _SUNDOMEIGEST_POWER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClassSUNDomEigEstimatorContent_Power_ =
  nb::class_<SUNDomEigEstimatorContent_Power_>(m, "SUNDomEigEstimatorContent_Power_",
                                               "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNDomEigEstimator_Power",
  [](N_Vector q, long max_iters, sunrealtype rel_tol, SUNContext sunctx)
    -> std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>>
  {
    auto SUNDomEigEstimator_Power_adapt_return_type_to_shared_ptr =
      [](N_Vector q, long max_iters, sunrealtype rel_tol, SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>>
    {
      auto lambda_result = SUNDomEigEstimator_Power(q, max_iters, rel_tol,
                                                    sunctx);

      return our_make_shared<std::remove_pointer_t<SUNDomEigEstimator>,
                             SUNDomEigEstimatorDeleter>(lambda_result);
    };

    return SUNDomEigEstimator_Power_adapt_return_type_to_shared_ptr(q, max_iters,
                                                                    rel_tol,
                                                                    sunctx);
  },
  nb::arg("q"), nb::arg("max_iters"), nb::arg("rel_tol"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def("SUNDomEigEstimator_SetRHS_Power", SUNDomEigEstimator_SetRHS_Power,
      nb::arg("DEE"), nb::arg("rhs_data"), nb::arg("RHSfn"));

m.def("SUNDomEigEstimator_SetRHSLinearizationVector_Power",
      SUNDomEigEstimator_SetRHSLinearizationVector_Power, nb::arg("DEE"),
      nb::arg("v"));

m.def("SUNDomEigEstimator_SetDEEisReal_Power",
      SUNDomEigEstimator_SetDEEisReal_Power, nb::arg("DEE"), nb::arg("real"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
