// #ifndef _SUNDOMEIGEST_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClassSUNDomEigEstimator_Ops_ =
  nb::class_<SUNDomEigEstimator_Ops_>(m,
                                      "SUNDomEigEstimator_Ops_", "Structure containing function pointers to estimator operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClassSUNDomEigEstimator_ =
  nb::class_<SUNDomEigEstimator_>(m,
                                  "SUNDomEigEstimator_", " An estimator is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of estimator\n   operations corresponding to that implementation.")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNDomEigEstimator_SetMaxIters", SUNDomEigEstimator_SetMaxIters,
      nb::arg("DEE"), nb::arg("max_iters"));

m.def("SUNDomEigEstimator_SetNumPreprocessIters",
      SUNDomEigEstimator_SetNumPreprocessIters, nb::arg("DEE"),
      nb::arg("num_iters"));

m.def("SUNDomEigEstimator_SetRelTol", SUNDomEigEstimator_SetRelTol,
      nb::arg("DEE"), nb::arg("tol"));

m.def("SUNDomEigEstimator_SetInitialGuess", SUNDomEigEstimator_SetInitialGuess,
      nb::arg("DEE"), nb::arg("q"));

m.def("SUNDomEigEstimator_Initialize", SUNDomEigEstimator_Initialize,
      nb::arg("DEE"));

m.def(
  "SUNDomEigEstimator_Estimate",
  [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, sunrealtype, sunrealtype>
  {
    auto SUNDomEigEstimator_Estimate_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, sunrealtype, sunrealtype>
    {
      sunrealtype lambdaR_adapt_modifiable;
      sunrealtype lambdaI_adapt_modifiable;

      SUNErrCode r = SUNDomEigEstimator_Estimate(DEE, &lambdaR_adapt_modifiable,
                                                 &lambdaI_adapt_modifiable);
      return std::make_tuple(r, lambdaR_adapt_modifiable,
                             lambdaI_adapt_modifiable);
    };

    return SUNDomEigEstimator_Estimate_adapt_modifiable_immutable_to_return(DEE);
  },
  nb::arg("DEE"));

m.def(
  "SUNDomEigEstimator_GetRes",
  [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNDomEigEstimator_GetRes_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype res_adapt_modifiable;

      SUNErrCode r = SUNDomEigEstimator_GetRes(DEE, &res_adapt_modifiable);
      return std::make_tuple(r, res_adapt_modifiable);
    };

    return SUNDomEigEstimator_GetRes_adapt_modifiable_immutable_to_return(DEE);
  },
  nb::arg("DEE"));

m.def(
  "SUNDomEigEstimator_GetNumIters",
  [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, long>
  {
    auto SUNDomEigEstimator_GetNumIters_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, long>
    {
      long num_iters_adapt_modifiable;

      SUNErrCode r =
        SUNDomEigEstimator_GetNumIters(DEE, &num_iters_adapt_modifiable);
      return std::make_tuple(r, num_iters_adapt_modifiable);
    };

    return SUNDomEigEstimator_GetNumIters_adapt_modifiable_immutable_to_return(
      DEE);
  },
  nb::arg("DEE"));

m.def(
  "SUNDomEigEstimator_GetNumATimesCalls",
  [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, long>
  {
    auto SUNDomEigEstimator_GetNumATimesCalls_adapt_modifiable_immutable_to_return =
      [](SUNDomEigEstimator DEE) -> std::tuple<SUNErrCode, long>
    {
      long num_ATimes_adapt_modifiable;

      SUNErrCode r =
        SUNDomEigEstimator_GetNumATimesCalls(DEE, &num_ATimes_adapt_modifiable);
      return std::make_tuple(r, num_ATimes_adapt_modifiable);
    };

    return SUNDomEigEstimator_GetNumATimesCalls_adapt_modifiable_immutable_to_return(
      DEE);
  },
  nb::arg("DEE"));

m.def("SUNDomEigEstimator_Write", SUNDomEigEstimator_Write, nb::arg("DEE"),
      nb::arg("outfile"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
