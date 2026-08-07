// #ifndef _SUNMATRIX_DENSE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_SUNMatrixContent_Dense =
  nb::class_<_SUNMatrixContent_Dense>(m, "_SUNMatrixContent_Dense", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNDenseMatrix",
  [](sunindextype M, sunindextype N,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNDenseMatrix_adapt_return_type_to_shared_ptr =
      [](sunindextype M, sunindextype N,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNDenseMatrix(M, N, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNDenseMatrix_adapt_return_type_to_shared_ptr(M, N, sunctx);
  },
  nb::arg("M"), nb::arg("N"), nb::arg("sunctx"), "nb::keep_alive<0, 3>()",
  nb::keep_alive<0, 3>());

m.def("SUNDenseMatrix_Print", SUNDenseMatrix_Print, nb::arg("A"),
      nb::arg("outfile"));

m.def("SUNDenseMatrix_Rows", SUNDenseMatrix_Rows, nb::arg("A"));

m.def("SUNDenseMatrix_Columns", SUNDenseMatrix_Columns, nb::arg("A"));

m.def("SUNDenseMatrix_LData", SUNDenseMatrix_LData, nb::arg("A"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
