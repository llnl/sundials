// #ifndef SUNDIALS_SUNDIALS_LOGGER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNLogLevel =
  nb::enum_<SUNLogLevel>(m, "SUNLogLevel", nb::is_arithmetic(), "")
    .value("SUN_LOGLEVEL_ALL", SUN_LOGLEVEL_ALL, "")
    .value("SUN_LOGLEVEL_NONE", SUN_LOGLEVEL_NONE, "")
    .value("SUN_LOGLEVEL_ERROR", SUN_LOGLEVEL_ERROR, "")
    .value("SUN_LOGLEVEL_WARNING", SUN_LOGLEVEL_WARNING, "")
    .value("SUN_LOGLEVEL_INFO", SUN_LOGLEVEL_INFO, "")
    .value("SUN_LOGLEVEL_DEBUG", SUN_LOGLEVEL_DEBUG, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def("SUNLogger_SetErrorFilename", SUNLogger_SetErrorFilename,
      nb::arg("logger"), nb::arg("error_filename"));

m.def("SUNLogger_SetWarningFilename", SUNLogger_SetWarningFilename,
      nb::arg("logger"), nb::arg("warning_filename"));

m.def("SUNLogger_SetDebugFilename", SUNLogger_SetDebugFilename,
      nb::arg("logger"), nb::arg("debug_filename"));

m.def("SUNLogger_SetInfoFilename", SUNLogger_SetInfoFilename, nb::arg("logger"),
      nb::arg("info_filename"));

m.def(
  "SUNLogger_QueueMsg",
  [](SUNLogger logger, SUNLogLevel lvl, const char* scope, const char* label,
     const char* msg_txt) -> SUNErrCode
  {
    auto SUNLogger_QueueMsg_adapt_variadic_format =
      [](SUNLogger logger, SUNLogLevel lvl, const char* scope,
         const char* label, const char* msg_txt) -> SUNErrCode
    {
      auto lambda_result = SUNLogger_QueueMsg(logger, lvl, scope, label, "%s",
                                              msg_txt);
      return lambda_result;
    };

    return SUNLogger_QueueMsg_adapt_variadic_format(logger, lvl, scope, label,
                                                    msg_txt);
  },
  nb::arg("logger"), nb::arg("lvl"), nb::arg("scope"), nb::arg("label"),
  nb::arg("msg_txt"));

m.def("SUNLogger_Flush", SUNLogger_Flush, nb::arg("logger"), nb::arg("lvl"));

m.def(
  "SUNLogger_GetOutputRank",
  [](SUNLogger logger) -> std::tuple<SUNErrCode, int>
  {
    auto SUNLogger_GetOutputRank_adapt_modifiable_immutable_to_return =
      [](SUNLogger logger) -> std::tuple<SUNErrCode, int>
    {
      int output_rank_adapt_modifiable;

      SUNErrCode r = SUNLogger_GetOutputRank(logger,
                                             &output_rank_adapt_modifiable);
      return std::make_tuple(r, output_rank_adapt_modifiable);
    };

    return SUNLogger_GetOutputRank_adapt_modifiable_immutable_to_return(logger);
  },
  nb::arg("logger"));
// #ifdef __cplusplus
// #endif
//
// #endif
