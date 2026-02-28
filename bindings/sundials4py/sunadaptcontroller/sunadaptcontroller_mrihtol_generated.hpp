// #ifndef _SUNADAPTCONTROLLER_MRIHTOL_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClassSUNAdaptControllerContent_MRIHTol_ =
  nb::class_<SUNAdaptControllerContent_MRIHTol_>(m, "SUNAdaptControllerContent_MRIHTol_",
                                                 "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNAdaptController_MRIHTol",
  [](SUNAdaptController HControl, SUNAdaptController TolControl,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_MRIHTol_adapt_return_type_to_shared_ptr =
      [](SUNAdaptController HControl, SUNAdaptController TolControl,
         SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_MRIHTol(HControl, TolControl,
                                                      sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_MRIHTol_adapt_return_type_to_shared_ptr(HControl,
                                                                      TolControl,
                                                                      sunctx);
  },
  nb::arg("HControl"), nb::arg("TolControl"), nb::arg("sunctx"),
  "nb::keep_alive<0, 3>()", nb::keep_alive<0, 3>());

m.def("SUNAdaptController_SetParams_MRIHTol",
      SUNAdaptController_SetParams_MRIHTol, nb::arg("C"),
      nb::arg("inner_max_relch"), nb::arg("inner_min_tolfac"),
      nb::arg("inner_max_tolfac"));

m.def(
  "SUNAdaptController_GetSlowController_MRIHTol",
  [](SUNAdaptController C) -> std::tuple<SUNErrCode, SUNAdaptController>
  {
    auto SUNAdaptController_GetSlowController_MRIHTol_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C) -> std::tuple<SUNErrCode, SUNAdaptController>
    {
      SUNAdaptController Cslow_adapt_modifiable;

      SUNErrCode r =
        SUNAdaptController_GetSlowController_MRIHTol(C, &Cslow_adapt_modifiable);
      return std::make_tuple(r, Cslow_adapt_modifiable);
    };

    return SUNAdaptController_GetSlowController_MRIHTol_adapt_modifiable_immutable_to_return(
      C);
  },
  nb::arg("C"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "SUNAdaptController_GetFastController_MRIHTol",
  [](SUNAdaptController C) -> std::tuple<SUNErrCode, SUNAdaptController>
  {
    auto SUNAdaptController_GetFastController_MRIHTol_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C) -> std::tuple<SUNErrCode, SUNAdaptController>
    {
      SUNAdaptController Cfast_adapt_modifiable;

      SUNErrCode r =
        SUNAdaptController_GetFastController_MRIHTol(C, &Cfast_adapt_modifiable);
      return std::make_tuple(r, Cfast_adapt_modifiable);
    };

    return SUNAdaptController_GetFastController_MRIHTol_adapt_modifiable_immutable_to_return(
      C);
  },
  nb::arg("C"), "nb::rv_policy::reference", nb::rv_policy::reference);
// #ifdef __cplusplus
//
// #endif
//
// #endif
