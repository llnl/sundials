// #ifndef _SUNMATRIX_BAND_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNMatrixContent_Band =
  nb::class_<_SUNMatrixContent_Band>(m, "_SUNMatrixContent_Band", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNBandMatrix",
  [](sunindextype N, sunindextype mu, sunindextype ml,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNBandMatrix_adapt_return_type_to_shared_ptr =
      [](sunindextype N, sunindextype mu, sunindextype ml,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNBandMatrix(N, mu, ml, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNBandMatrix_adapt_return_type_to_shared_ptr(N, mu, ml, sunctx);
  },
  nb::arg("N"), nb::arg("mu"), nb::arg("ml"), nb::arg("sunctx"),
  "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def(
  "SUNBandMatrixStorage",
  [](sunindextype N, sunindextype mu, sunindextype ml, sunindextype smu,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNBandMatrixStorage_adapt_return_type_to_shared_ptr =
      [](sunindextype N, sunindextype mu, sunindextype ml, sunindextype smu,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNBandMatrixStorage(N, mu, ml, smu, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNBandMatrixStorage_adapt_return_type_to_shared_ptr(N, mu, ml, smu,
                                                                sunctx);
  },
  nb::arg("N"), nb::arg("mu"), nb::arg("ml"), nb::arg("smu"), nb::arg("sunctx"),
  "nb::keep_alive<0, 5>()", nb::keep_alive<0, 5>());

m.def("SUNBandMatrix_Print", SUNBandMatrix_Print, nb::arg("A"),
      nb::arg("outfile"));

m.def("SUNBandMatrix_Rows", SUNBandMatrix_Rows, nb::arg("A"));

m.def("SUNBandMatrix_Columns", SUNBandMatrix_Columns, nb::arg("A"));

m.def("SUNBandMatrix_LowerBandwidth", SUNBandMatrix_LowerBandwidth, nb::arg("A"));

m.def("SUNBandMatrix_UpperBandwidth", SUNBandMatrix_UpperBandwidth, nb::arg("A"));

m.def("SUNBandMatrix_StoredUpperBandwidth", SUNBandMatrix_StoredUpperBandwidth,
      nb::arg("A"));

m.def("SUNBandMatrix_LDim", SUNBandMatrix_LDim, nb::arg("A"));

m.def("SUNBandMatrix_LData", SUNBandMatrix_LData, nb::arg("A"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
