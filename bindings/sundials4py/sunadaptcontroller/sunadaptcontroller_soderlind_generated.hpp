// #ifndef _SUNADAPTCONTROLLER_SODERLIND_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNAdaptControllerContent_Soderlind =
  nb::class_<_SUNAdaptControllerContent_Soderlind>(m, "_SUNAdaptControllerContent_Soderlind",
                                                   "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNAdaptController_Soderlind",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_Soderlind_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_Soderlind(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_Soderlind_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_Soderlind",
      SUNAdaptController_SetParams_Soderlind, nb::arg("C"), nb::arg("k1"),
      nb::arg("k2"), nb::arg("k3"), nb::arg("k4"), nb::arg("k5"));

m.def(
  "SUNAdaptController_PID",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_PID_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_PID(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_PID_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_PID", SUNAdaptController_SetParams_PID,
      nb::arg("C"), nb::arg("k1"), nb::arg("k2"), nb::arg("k3"));

m.def(
  "SUNAdaptController_PI",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_PI_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_PI(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_PI_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_PI", SUNAdaptController_SetParams_PI,
      nb::arg("C"), nb::arg("k1"), nb::arg("k2"));

m.def(
  "SUNAdaptController_I",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_I_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_I(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_I_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_I", SUNAdaptController_SetParams_I,
      nb::arg("C"), nb::arg("k1"));

m.def(
  "SUNAdaptController_ExpGus",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_ExpGus_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_ExpGus(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_ExpGus_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_ExpGus", SUNAdaptController_SetParams_ExpGus,
      nb::arg("C"), nb::arg("k1"), nb::arg("k2"));

m.def(
  "SUNAdaptController_ImpGus",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_ImpGus_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_ImpGus(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_ImpGus_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def("SUNAdaptController_SetParams_ImpGus", SUNAdaptController_SetParams_ImpGus,
      nb::arg("C"), nb::arg("k1"), nb::arg("k2"));

m.def(
  "SUNAdaptController_H0211",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_H0211_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_H0211(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_H0211_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def(
  "SUNAdaptController_H0321",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_H0321_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_H0321(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_H0321_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def(
  "SUNAdaptController_H211",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_H211_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_H211(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_H211_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def(
  "SUNAdaptController_H312",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
  {
    auto SUNAdaptController_H312_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNAdaptController>>
    {
      auto lambda_result = SUNAdaptController_H312(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNAdaptController>,
                             SUNAdaptControllerDeleter>(lambda_result);
    };

    return SUNAdaptController_H312_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
