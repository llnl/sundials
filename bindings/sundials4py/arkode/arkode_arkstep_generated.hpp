// #ifndef _ARKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

m.def("ARKStepSetExplicit", ARKStepSetExplicit, nb::arg("arkode_mem"));

m.def("ARKStepSetImplicit", ARKStepSetImplicit, nb::arg("arkode_mem"));

m.def("ARKStepSetImEx", ARKStepSetImEx, nb::arg("arkode_mem"));

m.def("ARKStepSetTables", ARKStepSetTables, nb::arg("arkode_mem"), nb::arg("q"),
      nb::arg("p"), nb::arg("Bi"), nb::arg("Be"));

m.def("ARKStepSetTableNum", ARKStepSetTableNum, nb::arg("arkode_mem"),
      nb::arg("itable"), nb::arg("etable"));

m.def("ARKStepSetTableName", ARKStepSetTableName, nb::arg("arkode_mem"),
      nb::arg("itable"), nb::arg("etable"));

m.def(
  "ARKStepGetCurrentButcherTables",
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

m.def(
  "ARKStepGetTimestepperStats",
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
