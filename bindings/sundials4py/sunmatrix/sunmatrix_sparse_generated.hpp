// #ifndef _SUNMATRIX_SPARSE_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("SUN_CSC_MAT") = 0;
m.attr("SUN_CSR_MAT") = 1;

auto pyClass_SUNMatrixContent_Sparse =
  nb::class_<_SUNMatrixContent_Sparse>(m, "_SUNMatrixContent_Sparse", "")
    .def(nb::init<>()) // implicit default constructor
  ;

sundials4py::scoped_def(
  m, "SUNSparseMatrix",
  [](sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNSparseMatrix_adapt_return_type_to_shared_ptr =
      [](sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNSparseMatrix(M, N, NNZ, sparsetype, sunctx);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNSparseMatrix_adapt_return_type_to_shared_ptr(M, N, NNZ,
                                                           sparsetype, sunctx);
  },
  nb::arg("M"), nb::arg("N"), nb::arg("NNZ"), nb::arg("sparsetype"),
  nb::arg("sunctx"), "nb::keep_alive<0, 5>()", nb::keep_alive<0, 5>());

sundials4py::scoped_def(
  m, "SUNSparseFromDenseMatrix",
  [](SUNMatrix A, sunrealtype droptol,
     int sparsetype) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNSparseFromDenseMatrix_adapt_return_type_to_shared_ptr =
      [](SUNMatrix A, sunrealtype droptol,
         int sparsetype) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNSparseFromDenseMatrix(A, droptol, sparsetype);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNSparseFromDenseMatrix_adapt_return_type_to_shared_ptr(A, droptol,
                                                                    sparsetype);
  },
  nb::arg("A"), nb::arg("droptol"), nb::arg("sparsetype"));

sundials4py::scoped_def(
  m, "SUNSparseFromBandMatrix",
  [](SUNMatrix A, sunrealtype droptol,
     int sparsetype) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNSparseFromBandMatrix_adapt_return_type_to_shared_ptr =
      [](SUNMatrix A, sunrealtype droptol,
         int sparsetype) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNSparseFromBandMatrix(A, droptol, sparsetype);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNSparseFromBandMatrix_adapt_return_type_to_shared_ptr(A, droptol,
                                                                   sparsetype);
  },
  nb::arg("A"), nb::arg("droptol"), nb::arg("sparsetype"));

sundials4py::scoped_def(
  m, "SUNSparseMatrix_ToCSR",
  [](const SUNMatrix A)
    -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNMatrix>>>
  {
    auto SUNSparseMatrix_ToCSR_adapt_modifiable_immutable_to_return =
      [](const SUNMatrix A) -> std::tuple<SUNErrCode, SUNMatrix>
    {
      SUNMatrix Bout_adapt_modifiable;

      SUNErrCode r = SUNSparseMatrix_ToCSR(A, &Bout_adapt_modifiable);
      return std::make_tuple(r, Bout_adapt_modifiable);
    };
    auto SUNSparseMatrix_ToCSR_adapt_return_type_to_shared_ptr =
      [&SUNSparseMatrix_ToCSR_adapt_modifiable_immutable_to_return](
        const SUNMatrix A)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNMatrix>>>
    {
      auto lambda_result =
        SUNSparseMatrix_ToCSR_adapt_modifiable_immutable_to_return(A);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNMatrix>,
                                             SUNMatrixDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SUNSparseMatrix_ToCSR_adapt_return_type_to_shared_ptr(A);
  },
  nb::arg("A"), nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "SUNSparseMatrix_ToCSC",
  [](const SUNMatrix A)
    -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNMatrix>>>
  {
    auto SUNSparseMatrix_ToCSC_adapt_modifiable_immutable_to_return =
      [](const SUNMatrix A) -> std::tuple<SUNErrCode, SUNMatrix>
    {
      SUNMatrix Bout_adapt_modifiable;

      SUNErrCode r = SUNSparseMatrix_ToCSC(A, &Bout_adapt_modifiable);
      return std::make_tuple(r, Bout_adapt_modifiable);
    };
    auto SUNSparseMatrix_ToCSC_adapt_return_type_to_shared_ptr =
      [&SUNSparseMatrix_ToCSC_adapt_modifiable_immutable_to_return](
        const SUNMatrix A)
      -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNMatrix>>>
    {
      auto lambda_result =
        SUNSparseMatrix_ToCSC_adapt_modifiable_immutable_to_return(A);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNMatrix>,
                                             SUNMatrixDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SUNSparseMatrix_ToCSC_adapt_return_type_to_shared_ptr(A);
  },
  nb::arg("A"), nb::rv_policy::reference);

sundials4py::scoped_def(m, "SUNSparseMatrix_Realloc", SUNSparseMatrix_Realloc,
                        nb::arg("A"));

sundials4py::scoped_def(m, "SUNSparseMatrix_Reallocate",
                        SUNSparseMatrix_Reallocate, nb::arg("A"), nb::arg("NNZ"));

sundials4py::scoped_def(m, "SUNSparseMatrix_Print", SUNSparseMatrix_Print,
                        nb::arg("A"), nb::arg("outfile"));

sundials4py::scoped_def(m, "SUNSparseMatrix_Rows", SUNSparseMatrix_Rows,
                        nb::arg("A"));

sundials4py::scoped_def(m, "SUNSparseMatrix_Columns", SUNSparseMatrix_Columns,
                        nb::arg("A"));

sundials4py::scoped_def(m, "SUNSparseMatrix_NNZ", SUNSparseMatrix_NNZ,
                        nb::arg("A"));

sundials4py::scoped_def(m, "SUNSparseMatrix_NP", SUNSparseMatrix_NP,
                        nb::arg("A"));

sundials4py::scoped_def(m, "SUNSparseMatrix_SparseType",
                        SUNSparseMatrix_SparseType, nb::arg("A"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
