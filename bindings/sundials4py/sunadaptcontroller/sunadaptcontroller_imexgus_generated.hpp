// #ifndef _SUNADAPTCONTROLLER_IMEXGUS_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNAdaptControllerContent_ImExGus =
  nb::class_<_SUNAdaptControllerContent_ImExGus>(m, "_SUNAdaptControllerContent_ImExGus",
                                                 "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNAdaptController_ImExGus",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_ImExGus_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_ImExGus(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_ImExGus_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_ImExGus",
      SUNAdaptController_SetParams_ImExGus, nb::arg("C"), nb::arg("k1e"),
      nb::arg("k2e"), nb::arg("k1i"), nb::arg("k2i"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
