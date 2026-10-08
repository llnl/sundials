// #ifndef _ARKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

sundials4py::scoped_def(m, "ARKStepSetExplicit", ARKStepSetExplicit,
                        nb::arg("arkode_mem"));

sundials4py::scoped_def(m, "ARKStepSetImplicit", ARKStepSetImplicit,
                        nb::arg("arkode_mem"));

sundials4py::scoped_def(m, "ARKStepSetImEx", ARKStepSetImEx,
                        nb::arg("arkode_mem"));

sundials4py::scoped_def(
  m, "ARKStepSetTables",
  [](void* arkode_mem, int q, int p,
     std::optional<ARKodeButcherTable> Bi = std::nullopt,
     std::optional<ARKodeButcherTable> Be = std::nullopt) -> int
  {
    auto ARKStepSetTables_adapt_optional_arg_with_default_null =
      [](void* arkode_mem, int q, int p,
         std::optional<ARKodeButcherTable> Bi = std::nullopt,
         std::optional<ARKodeButcherTable> Be = std::nullopt) -> int
    {
      ARKodeButcherTable Bi_adapt_default_null = nullptr;
      if (Bi.has_value()) Bi_adapt_default_null = Bi.value();
      ARKodeButcherTable Be_adapt_default_null = nullptr;
      if (Be.has_value()) Be_adapt_default_null = Be.value();

      auto lambda_result = ARKStepSetTables(arkode_mem, q, p,
                                            Bi_adapt_default_null,
                                            Be_adapt_default_null);
      return lambda_result;
    };

    return ARKStepSetTables_adapt_optional_arg_with_default_null(arkode_mem, q,
                                                                 p, Bi, Be);
  },
  nb::arg("arkode_mem"), nb::arg("q"), nb::arg("p"),
  nb::arg("Bi").none() = nb::none(), nb::arg("Be").none() = nb::none());

sundials4py::scoped_def(m, "ARKStepSetTableNum", ARKStepSetTableNum,
                        nb::arg("arkode_mem"), nb::arg("itable"),
                        nb::arg("etable"));

sundials4py::scoped_def(m, "ARKStepSetTableName", ARKStepSetTableName,
                        nb::arg("arkode_mem"), nb::arg("itable"),
                        nb::arg("etable"));

sundials4py::scoped_def(
  m, "ARKStepGetCurrentButcherTables",
  [](void* arkode_mem) -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
  {
    auto ARKStepGetCurrentButcherTables_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
    {
      ARKodeButcherTable Bi_adapt_modifiable;
      ARKodeButcherTable Be_adapt_modifiable;

      int r = ARKStepGetCurrentButcherTables(arkode_mem, &Bi_adapt_modifiable,
                                             &Be_adapt_modifiable);
      return std::make_tuple(r, Bi_adapt_modifiable, Be_adapt_modifiable);
    };

    return ARKStepGetCurrentButcherTables_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"),
  " Optional output functions\n\n nb::rv_policy::reference",
  nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ARKStepGetTimestepperStats",
  [](void* arkode_mem) -> std::tuple<int, long, long, long, long, long, long, long>
  {
    auto ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem)
      -> std::tuple<int, long, long, long, long, long, long, long>
    {
      long expsteps_adapt_modifiable;
      long accsteps_adapt_modifiable;
      long step_attempts_adapt_modifiable;
      long nfe_evals_adapt_modifiable;
      long nfi_evals_adapt_modifiable;
      long nlinsetups_adapt_modifiable;
      long netfails_adapt_modifiable;

      int r = ARKStepGetTimestepperStats(arkode_mem, &expsteps_adapt_modifiable,
                                         &accsteps_adapt_modifiable,
                                         &step_attempts_adapt_modifiable,
                                         &nfe_evals_adapt_modifiable,
                                         &nfi_evals_adapt_modifiable,
                                         &nlinsetups_adapt_modifiable,
                                         &netfails_adapt_modifiable);
      return std::make_tuple(r, expsteps_adapt_modifiable,
                             accsteps_adapt_modifiable,
                             step_attempts_adapt_modifiable,
                             nfe_evals_adapt_modifiable,
                             nfi_evals_adapt_modifiable,
                             nlinsetups_adapt_modifiable,
                             netfails_adapt_modifiable);
    };

    return ARKStepGetTimestepperStats_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
