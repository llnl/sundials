// #ifndef _SUNMATRIX_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNMatrix_ID =
  nb::enum_<SUNMatrix_ID>(m, "SUNMatrix_ID", nb::is_arithmetic(), "")
    .value("SUNMATRIX_DENSE", SUNMATRIX_DENSE, "")
    .value("SUNMATRIX_MAGMADENSE", SUNMATRIX_MAGMADENSE, "")
    .value("SUNMATRIX_ONEMKLDENSE", SUNMATRIX_ONEMKLDENSE, "")
    .value("SUNMATRIX_BAND", SUNMATRIX_BAND, "")
    .value("SUNMATRIX_SPARSE", SUNMATRIX_SPARSE, "")
    .value("SUNMATRIX_SLUNRLOC", SUNMATRIX_SLUNRLOC, "")
    .value("SUNMATRIX_CUSPARSE", SUNMATRIX_CUSPARSE, "")
    .value("SUNMATRIX_GINKGO", SUNMATRIX_GINKGO, "")
    .value("SUNMATRIX_GINKGOBATCH", SUNMATRIX_GINKGOBATCH, "")
    .value("SUNMATRIX_KOKKOSDENSE", SUNMATRIX_KOKKOSDENSE, "")
    .value("SUNMATRIX_CUSTOM", SUNMATRIX_CUSTOM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNMatrix_Ops =
  nb::class_<_generic_SUNMatrix_Ops>(m,
                                     "_generic_SUNMatrix_Ops", "Structure containing function pointers to matrix operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNMatrix =
  nb::class_<_generic_SUNMatrix>(m,
                                 "_generic_SUNMatrix", " A matrix is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of matrix\n   operations corresponding to that implementation.")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNMatGetID", SUNMatGetID, nb::arg("A"));

m.def(
  "SUNMatClone",
  [](SUNMatrix A) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
  {
    auto SUNMatClone_adapt_return_type_to_shared_ptr =
      [](SUNMatrix A) -> std::shared_ptr<std::remove_pointer_t<SUNMatrix>>
    {
      auto lambda_result = SUNMatClone(A);

      return our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
        lambda_result);
    };

    return SUNMatClone_adapt_return_type_to_shared_ptr(A);
  },
  nb::arg("A"));

m.def("SUNMatZero", SUNMatZero, nb::arg("A"));

m.def("SUNMatCopy", SUNMatCopy, nb::arg("A"), nb::arg("B"));

m.def("SUNMatScaleAdd", SUNMatScaleAdd, nb::arg("c"), nb::arg("A"), nb::arg("B"));

m.def("SUNMatScaleAddI", SUNMatScaleAddI, nb::arg("c"), nb::arg("A"));

m.def("SUNMatMatvecSetup", SUNMatMatvecSetup, nb::arg("A"));

m.def("SUNMatMatvec", SUNMatMatvec, nb::arg("A"), nb::arg("x"), nb::arg("y"));

m.def("SUNMatHermitianTransposeVec", SUNMatHermitianTransposeVec, nb::arg("A"),
      nb::arg("x"), nb::arg("y"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
