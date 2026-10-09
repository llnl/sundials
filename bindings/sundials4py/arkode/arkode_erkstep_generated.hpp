// #ifndef _ERKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

sundials4py::scoped_def(m, "ERKStepSetTable", ERKStepSetTable,
                        nb::arg("arkode_mem"), nb::arg("B"));

sundials4py::scoped_def(m, "ERKStepSetTableNum", ERKStepSetTableNum,
                        nb::arg("arkode_mem"), nb::arg("etable"));

sundials4py::scoped_def(m, "ERKStepSetTableName", ERKStepSetTableName,
                        nb::arg("arkode_mem"), nb::arg("etable"));

sundials4py::scoped_def(
  m, "ERKStepGetCurrentButcherTable",
  [](void* arkode_mem) -> std::tuple<int, ARKodeButcherTable>
  {
    auto ERKStepGetCurrentButcherTable_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, ARKodeButcherTable>
    {
      ARKodeButcherTable B_adapt_modifiable;

      int r = ERKStepGetCurrentButcherTable(arkode_mem, &B_adapt_modifiable);
      return std::make_tuple(r, B_adapt_modifiable);
    };

    return ERKStepGetCurrentButcherTable_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"),
  " Optional output functions\n\n nb::rv_policy::reference",
  nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "ERKStepGetTimestepperStats",
  [](void* arkode_mem) -> std::tuple<int, long, long, long, long, long>
  {
    auto ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long, long, long, long, long>
    {
      long expsteps_adapt_modifiable;
      long accsteps_adapt_modifiable;
      long step_attempts_adapt_modifiable;
      long nfevals_adapt_modifiable;
      long netfails_adapt_modifiable;

      int r = ERKStepGetTimestepperStats(arkode_mem, &expsteps_adapt_modifiable,
                                         &accsteps_adapt_modifiable,
                                         &step_attempts_adapt_modifiable,
                                         &nfevals_adapt_modifiable,
                                         &netfails_adapt_modifiable);
      return std::make_tuple(r, expsteps_adapt_modifiable,
                             accsteps_adapt_modifiable,
                             step_attempts_adapt_modifiable,
                             nfevals_adapt_modifiable, netfails_adapt_modifiable);
    };

    return ERKStepGetTimestepperStats_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"), "Grouped optional output functions");
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
