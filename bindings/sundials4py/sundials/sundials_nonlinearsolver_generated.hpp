// #ifndef _SUNNONLINEARSOLVER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNNonlinearSolver_Type =
  nb::enum_<SUNNonlinearSolver_Type>(m, "SUNNonlinearSolver_Type",
                                     nb::is_arithmetic(), "")
    .value("SUNNONLINEARSOLVER_ROOTFIND", SUNNONLINEARSOLVER_ROOTFIND, "")
    .value("SUNNONLINEARSOLVER_FIXEDPOINT", SUNNONLINEARSOLVER_FIXEDPOINT, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNNonlinearSolver_Ops =
  nb::class_<
    _generic_SUNNonlinearSolver_Ops>(m, "_generic_SUNNonlinearSolver_Ops",
                                     "Structure containing function pointers "
                                     "to nonlinear solver operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNNonlinearSolver =
  nb::class_<_generic_SUNNonlinearSolver>(m,
                                          "_generic_SUNNonlinearSolver", " A nonlinear solver is a structure with an implementation-dependent 'content'\n   field, and a pointer to a structure of solver nonlinear solver operations\n   corresponding to that implementation.")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNNonlinSolGetType", SUNNonlinSolGetType, nb::arg("NLS"));

m.def("SUNNonlinSolInitialize", SUNNonlinSolInitialize, nb::arg("NLS"));

m.def("SUNNonlinSolSetMaxIters", SUNNonlinSolSetMaxIters, nb::arg("NLS"),
      nb::arg("maxiters"));

m.def(
  "SUNNonlinSolGetNumIters",
  [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long>
  {
    auto SUNNonlinSolGetNumIters_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long>
    {
      long niters_adapt_modifiable;

      SUNErrCode r = SUNNonlinSolGetNumIters(NLS, &niters_adapt_modifiable);
      return std::make_tuple(r, niters_adapt_modifiable);
    };

    return SUNNonlinSolGetNumIters_adapt_modifiable_immutable_to_return(NLS);
  },
  nb::arg("NLS"));

m.def(
  "SUNNonlinSolGetCurIter",
  [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, int>
  {
    auto SUNNonlinSolGetCurIter_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, int>
    {
      int iter_adapt_modifiable;

      SUNErrCode r = SUNNonlinSolGetCurIter(NLS, &iter_adapt_modifiable);
      return std::make_tuple(r, iter_adapt_modifiable);
    };

    return SUNNonlinSolGetCurIter_adapt_modifiable_immutable_to_return(NLS);
  },
  nb::arg("NLS"));

m.def(
  "SUNNonlinSolGetNumConvFails",
  [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long>
  {
    auto SUNNonlinSolGetNumConvFails_adapt_modifiable_immutable_to_return =
      [](SUNNonlinearSolver NLS) -> std::tuple<SUNErrCode, long>
    {
      long nconvfails_adapt_modifiable;

      SUNErrCode r = SUNNonlinSolGetNumConvFails(NLS,
                                                 &nconvfails_adapt_modifiable);
      return std::make_tuple(r, nconvfails_adapt_modifiable);
    };

    return SUNNonlinSolGetNumConvFails_adapt_modifiable_immutable_to_return(NLS);
  },
  nb::arg("NLS"));
m.attr("SUN_NLS_CONTINUE")   = +901;
m.attr("SUN_NLS_CONV_RECVR") = +902;
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
