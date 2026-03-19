// #ifndef ARKODE_FORCINGINGSTEP_H_
//
// #ifdef __cplusplus
// #endif
//

m.def("ForcingStepReInit", ForcingStepReInit, nb::arg("arkode_mem"),
      nb::arg("stepper1"), nb::arg("stepper2"), nb::arg("t0"), nb::arg("y0"));

m.def(
  "ForcingStepGetNumEvolves",
  [](void* arkode_mem, int partition) -> std::tuple<int, long>
  {
    auto ForcingStepGetNumEvolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int partition) -> std::tuple<int, long>
    {
      long evolves_adapt_modifiable;

      int r = ForcingStepGetNumEvolves(arkode_mem, partition,
                                       &evolves_adapt_modifiable);
      return std::make_tuple(r, evolves_adapt_modifiable);
    };

    return ForcingStepGetNumEvolves_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                         partition);
  },
  nb::arg("arkode_mem"), nb::arg("partition"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
