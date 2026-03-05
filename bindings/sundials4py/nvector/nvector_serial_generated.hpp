// #ifndef _NVECTOR_SERIAL_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_Serial =
  nb::class_<_N_VectorContent_Serial>(m, "_N_VectorContent_Serial", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "N_VNewEmpty_Serial",
  [](sunindextype vec_length,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VNewEmpty_Serial_adapt_return_type_to_shared_ptr =
      [](sunindextype vec_length,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VNewEmpty_Serial(vec_length, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VNewEmpty_Serial_adapt_return_type_to_shared_ptr(vec_length, sunctx);
  },
  nb::arg("vec_length"), nb::arg("sunctx"), "nb::keep_alive<0, 2>()",
  nb::keep_alive<0, 2>());

m.def(
  "N_VNew_Serial",
  [](sunindextype vec_length,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VNew_Serial_adapt_return_type_to_shared_ptr =
      [](sunindextype vec_length,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VNew_Serial(vec_length, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VNew_Serial_adapt_return_type_to_shared_ptr(vec_length, sunctx);
  },
  nb::arg("vec_length"), nb::arg("sunctx"), "nb::keep_alive<0, 2>()",
  nb::keep_alive<0, 2>());

m.def(
  "N_VMake_Serial",
  [](sunindextype vec_length, sundials4py::Array1d v_data_1d,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VMake_Serial_adapt_arr_ptr_to_std_vector =
      [](sunindextype vec_length, sundials4py::Array1d v_data_1d,
         SUNContext sunctx) -> N_Vector
    {
      sunrealtype* v_data_1d_ptr =
        reinterpret_cast<sunrealtype*>(v_data_1d.data());

      auto lambda_result = N_VMake_Serial(vec_length, v_data_1d_ptr, sunctx);
      return lambda_result;
    };
    auto N_VMake_Serial_adapt_return_type_to_shared_ptr =
      [&N_VMake_Serial_adapt_arr_ptr_to_std_vector](sunindextype vec_length,
                                                    sundials4py::Array1d v_data_1d,
                                                    SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result =
        N_VMake_Serial_adapt_arr_ptr_to_std_vector(vec_length, v_data_1d, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VMake_Serial_adapt_return_type_to_shared_ptr(vec_length, v_data_1d,
                                                          sunctx);
  },
  nb::arg("vec_length"), nb::arg("v_data_1d"), nb::arg("sunctx"),
  "nb::keep_alive<0, 3>()", nb::keep_alive<0, 3>());

m.def("N_VEnableFusedOps_Serial", N_VEnableFusedOps_Serial, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearCombination_Serial", N_VEnableLinearCombination_Serial,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMulti_Serial", N_VEnableScaleAddMulti_Serial,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableDotProdMulti_Serial", N_VEnableDotProdMulti_Serial,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearSumVectorArray_Serial",
      N_VEnableLinearSumVectorArray_Serial, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleVectorArray_Serial", N_VEnableScaleVectorArray_Serial,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableConstVectorArray_Serial", N_VEnableConstVectorArray_Serial,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormVectorArray_Serial",
      N_VEnableWrmsNormVectorArray_Serial, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormMaskVectorArray_Serial",
      N_VEnableWrmsNormMaskVectorArray_Serial, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMultiVectorArray_Serial",
      N_VEnableScaleAddMultiVectorArray_Serial, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearCombinationVectorArray_Serial",
      N_VEnableLinearCombinationVectorArray_Serial, nb::arg("v"), nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
