// #ifndef _SUNDIALS_LOGGER_H
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

m.def(
  "SUNLogger_Create",
  [](SUNComm comm, int output_rank)
    -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNLogger>>>
  {
    auto SUNLogger_Create_adapt_modifiable_immutable_to_return =
      [](SUNComm comm, int output_rank) -> std::tuple<SUNErrCode, SUNLogger>
    {
      SUNLogger logger_adapt_modifiable;

      SUNErrCode r = SUNLogger_Create(comm, output_rank,
                                      &logger_adapt_modifiable);
      return std::make_tuple(r, logger_adapt_modifiable);
    };
    auto SUNLogger_Create_adapt_return_type_to_shared_ptr =
      [&SUNLogger_Create_adapt_modifiable_immutable_to_return](SUNComm comm,
                                                               int output_rank)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNLogger>>>
    {
      auto lambda_result =
        SUNLogger_Create_adapt_modifiable_immutable_to_return(comm, output_rank);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNLogger>,
                                             SUNLoggerDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SUNLogger_Create_adapt_return_type_to_shared_ptr(comm, output_rank);
  },
  nb::arg("comm"), nb::arg("output_rank"), nb::rv_policy::reference);

m.def(
  "SUNLogger_CreateFromEnv",
  [](SUNComm comm)
    -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNLogger>>>
  {
    auto SUNLogger_CreateFromEnv_adapt_modifiable_immutable_to_return =
      [](SUNComm comm) -> std::tuple<SUNErrCode, SUNLogger>
    {
      SUNLogger logger_adapt_modifiable;

      SUNErrCode r = SUNLogger_CreateFromEnv(comm, &logger_adapt_modifiable);
      return std::make_tuple(r, logger_adapt_modifiable);
    };
    auto SUNLogger_CreateFromEnv_adapt_return_type_to_shared_ptr =
      [&SUNLogger_CreateFromEnv_adapt_modifiable_immutable_to_return](SUNComm comm)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNLogger>>>
    {
      auto lambda_result =
        SUNLogger_CreateFromEnv_adapt_modifiable_immutable_to_return(comm);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNLogger>,
                                             SUNLoggerDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SUNLogger_CreateFromEnv_adapt_return_type_to_shared_ptr(comm);
  },
  nb::arg("comm"), nb::rv_policy::reference);

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
