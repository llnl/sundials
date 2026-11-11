// #ifndef SUNDIALS_SUNDIALS_CONTEXT_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SUNContext_GetLastError", SUNContext_GetLastError, nb::arg("sunctx"));

m.def("SUNContext_PeekLastError", SUNContext_PeekLastError, nb::arg("sunctx"));

m.def("SUNContext_ClearErrHandlers", SUNContext_ClearErrHandlers,
      nb::arg("sunctx"));

m.def(
  "SUNContext_GetProfiler",
  [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNProfiler>
  {
    auto SUNContext_GetProfiler_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNProfiler>
    {
      SUNProfiler profiler_adapt_modifiable;

      SUNErrCode r = SUNContext_GetProfiler(sunctx, &profiler_adapt_modifiable);
      return std::make_tuple(r, profiler_adapt_modifiable);
    };

    return SUNContext_GetProfiler_adapt_modifiable_immutable_to_return(sunctx);
  },
  nb::arg("sunctx"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def("SUNContext_SetProfiler", SUNContext_SetProfiler, nb::arg("sunctx"),
      nb::arg("profiler"));

m.def(
  "SUNContext_GetLogger",
  [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNLogger>
  {
    auto SUNContext_GetLogger_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNLogger>
    {
      SUNLogger logger_adapt_modifiable;

      SUNErrCode r = SUNContext_GetLogger(sunctx, &logger_adapt_modifiable);
      return std::make_tuple(r, logger_adapt_modifiable);
    };

    return SUNContext_GetLogger_adapt_modifiable_immutable_to_return(sunctx);
  },
  nb::arg("sunctx"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def("SUNContext_SetLogger", SUNContext_SetLogger, nb::arg("sunctx"),
      nb::arg("logger"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
