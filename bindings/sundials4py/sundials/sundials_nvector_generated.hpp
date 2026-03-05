// #ifndef _NVECTOR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumN_Vector_ID =
  nb::enum_<N_Vector_ID>(m, "N_Vector_ID", nb::is_arithmetic(), "")
    .value("SUNDIALS_NVEC_SERIAL", SUNDIALS_NVEC_SERIAL, "")
    .value("SUNDIALS_NVEC_PARALLEL", SUNDIALS_NVEC_PARALLEL, "")
    .value("SUNDIALS_NVEC_OPENMP", SUNDIALS_NVEC_OPENMP, "")
    .value("SUNDIALS_NVEC_PTHREADS", SUNDIALS_NVEC_PTHREADS, "")
    .value("SUNDIALS_NVEC_PARHYP", SUNDIALS_NVEC_PARHYP, "")
    .value("SUNDIALS_NVEC_PETSC", SUNDIALS_NVEC_PETSC, "")
    .value("SUNDIALS_NVEC_CUDA", SUNDIALS_NVEC_CUDA, "")
    .value("SUNDIALS_NVEC_HIP", SUNDIALS_NVEC_HIP, "")
    .value("SUNDIALS_NVEC_SYCL", SUNDIALS_NVEC_SYCL, "")
    .value("SUNDIALS_NVEC_RAJA", SUNDIALS_NVEC_RAJA, "")
    .value("SUNDIALS_NVEC_KOKKOS", SUNDIALS_NVEC_KOKKOS, "")
    .value("SUNDIALS_NVEC_OPENMPDEV", SUNDIALS_NVEC_OPENMPDEV, "")
    .value("SUNDIALS_NVEC_TRILINOS", SUNDIALS_NVEC_TRILINOS, "")
    .value("SUNDIALS_NVEC_MANYVECTOR", SUNDIALS_NVEC_MANYVECTOR, "")
    .value("SUNDIALS_NVEC_MPIMANYVECTOR", SUNDIALS_NVEC_MPIMANYVECTOR, "")
    .value("SUNDIALS_NVEC_MPIPLUSX", SUNDIALS_NVEC_MPIPLUSX, "")
    .value("SUNDIALS_NVEC_CUSTOM", SUNDIALS_NVEC_CUSTOM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_N_Vector_Ops =
  nb::class_<_generic_N_Vector_Ops>(m,
                                    "_generic_N_Vector_Ops", "Structure containing function pointers to vector operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_N_Vector =
  nb::class_<_generic_N_Vector>(m,
                                "_generic_N_Vector", " A vector is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of vector\n   operations corresponding to that implementation.")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("N_VGetVectorID", N_VGetVectorID, nb::arg("w"));

m.def(
  "N_VClone",
  [](N_Vector w) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VClone_adapt_return_type_to_shared_ptr =
      [](N_Vector w) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VClone(w);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VClone_adapt_return_type_to_shared_ptr(w);
  },
  nb::arg("w"));

m.def(
  "N_VCloneEmpty",
  [](N_Vector w) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VCloneEmpty_adapt_return_type_to_shared_ptr =
      [](N_Vector w) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VCloneEmpty(w);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VCloneEmpty_adapt_return_type_to_shared_ptr(w);
  },
  nb::arg("w"));

m.def("N_VGetCommunicator", N_VGetCommunicator, nb::arg("v"));

m.def("N_VGetLength", N_VGetLength, nb::arg("v"));

m.def("N_VGetLocalLength", N_VGetLocalLength, nb::arg("v"));

m.def("N_VLinearSum", N_VLinearSum, nb::arg("a"), nb::arg("x"), nb::arg("b"),
      nb::arg("y"), nb::arg("z"));

m.def("N_VConst", N_VConst, nb::arg("c"), nb::arg("z"));

m.def("N_VProd", N_VProd, nb::arg("x"), nb::arg("y"), nb::arg("z"));

m.def("N_VDiv", N_VDiv, nb::arg("x"), nb::arg("y"), nb::arg("z"));

m.def("N_VScale", N_VScale, nb::arg("c"), nb::arg("x"), nb::arg("z"));

m.def("N_VAbs", N_VAbs, nb::arg("x"), nb::arg("z"));

m.def("N_VInv", N_VInv, nb::arg("x"), nb::arg("z"));

m.def("N_VAddConst", N_VAddConst, nb::arg("x"), nb::arg("b"), nb::arg("z"));

m.def("N_VDotProd", N_VDotProd, nb::arg("x"), nb::arg("y"));

m.def("N_VMaxNorm", N_VMaxNorm, nb::arg("x"));

m.def("N_VWrmsNorm", N_VWrmsNorm, nb::arg("x"), nb::arg("w"));

m.def("N_VWrmsNormMask", N_VWrmsNormMask, nb::arg("x"), nb::arg("w"),
      nb::arg("id"));

m.def("N_VMin", N_VMin, nb::arg("x"));

m.def("N_VWL2Norm", N_VWL2Norm, nb::arg("x"), nb::arg("w"));

m.def("N_VL1Norm", N_VL1Norm, nb::arg("x"));

m.def("N_VCompare", N_VCompare, nb::arg("c"), nb::arg("x"), nb::arg("z"));

m.def("N_VInvTest", N_VInvTest, nb::arg("x"), nb::arg("z"));

m.def("N_VConstrMask", N_VConstrMask, nb::arg("c"), nb::arg("x"), nb::arg("m"));

m.def("N_VMinQuotient", N_VMinQuotient, nb::arg("num"), nb::arg("denom"));

m.def(
  "N_VLinearCombination",
  [](int nvec, sundials4py::Array1d c_1d, std::vector<N_Vector> X_1d,
     N_Vector z) -> SUNErrCode
  {
    auto N_VLinearCombination_adapt_arr_ptr_to_std_vector =
      [](int nvec, sundials4py::Array1d c_1d, std::vector<N_Vector> X_1d,
         N_Vector z) -> SUNErrCode
    {
      sunrealtype* c_1d_ptr = reinterpret_cast<sunrealtype*>(c_1d.data());
      N_Vector* X_1d_ptr =
        reinterpret_cast<N_Vector*>(X_1d.empty() ? nullptr : X_1d.data());

      auto lambda_result = N_VLinearCombination(nvec, c_1d_ptr, X_1d_ptr, z);
      return lambda_result;
    };

    return N_VLinearCombination_adapt_arr_ptr_to_std_vector(nvec, c_1d, X_1d, z);
  },
  nb::arg("nvec"), nb::arg("c_1d"), nb::arg("X_1d"), nb::arg("z"));

m.def(
  "N_VScaleAddMulti",
  [](int nvec, sundials4py::Array1d a_1d, N_Vector x,
     std::vector<N_Vector> Y_1d, std::vector<N_Vector> Z_1d) -> SUNErrCode
  {
    auto N_VScaleAddMulti_adapt_arr_ptr_to_std_vector =
      [](int nvec, sundials4py::Array1d a_1d, N_Vector x,
         std::vector<N_Vector> Y_1d, std::vector<N_Vector> Z_1d) -> SUNErrCode
    {
      sunrealtype* a_1d_ptr = reinterpret_cast<sunrealtype*>(a_1d.data());
      N_Vector* Y_1d_ptr =
        reinterpret_cast<N_Vector*>(Y_1d.empty() ? nullptr : Y_1d.data());
      N_Vector* Z_1d_ptr =
        reinterpret_cast<N_Vector*>(Z_1d.empty() ? nullptr : Z_1d.data());

      auto lambda_result = N_VScaleAddMulti(nvec, a_1d_ptr, x, Y_1d_ptr,
                                            Z_1d_ptr);
      return lambda_result;
    };

    return N_VScaleAddMulti_adapt_arr_ptr_to_std_vector(nvec, a_1d, x, Y_1d,
                                                        Z_1d);
  },
  nb::arg("nvec"), nb::arg("a_1d"), nb::arg("x"), nb::arg("Y_1d"),
  nb::arg("Z_1d"));

m.def(
  "N_VDotProdMulti",
  [](int nvec, N_Vector x, std::vector<N_Vector> Y_1d,
     sundials4py::Array1d dotprods_1d) -> SUNErrCode
  {
    auto N_VDotProdMulti_adapt_arr_ptr_to_std_vector =
      [](int nvec, N_Vector x, std::vector<N_Vector> Y_1d,
         sundials4py::Array1d dotprods_1d) -> SUNErrCode
    {
      N_Vector* Y_1d_ptr =
        reinterpret_cast<N_Vector*>(Y_1d.empty() ? nullptr : Y_1d.data());
      sunrealtype* dotprods_1d_ptr =
        reinterpret_cast<sunrealtype*>(dotprods_1d.data());

      auto lambda_result = N_VDotProdMulti(nvec, x, Y_1d_ptr, dotprods_1d_ptr);
      return lambda_result;
    };

    return N_VDotProdMulti_adapt_arr_ptr_to_std_vector(nvec, x, Y_1d,
                                                       dotprods_1d);
  },
  nb::arg("nvec"), nb::arg("x"), nb::arg("Y_1d"), nb::arg("dotprods_1d"));

m.def(
  "N_VLinearSumVectorArray",
  [](int nvec, sunrealtype a, std::vector<N_Vector> X_1d, sunrealtype b,
     std::vector<N_Vector> Y_1d, std::vector<N_Vector> Z_1d) -> SUNErrCode
  {
    auto N_VLinearSumVectorArray_adapt_arr_ptr_to_std_vector =
      [](int nvec, sunrealtype a, std::vector<N_Vector> X_1d, sunrealtype b,
         std::vector<N_Vector> Y_1d, std::vector<N_Vector> Z_1d) -> SUNErrCode
    {
      N_Vector* X_1d_ptr =
        reinterpret_cast<N_Vector*>(X_1d.empty() ? nullptr : X_1d.data());
      N_Vector* Y_1d_ptr =
        reinterpret_cast<N_Vector*>(Y_1d.empty() ? nullptr : Y_1d.data());
      N_Vector* Z_1d_ptr =
        reinterpret_cast<N_Vector*>(Z_1d.empty() ? nullptr : Z_1d.data());

      auto lambda_result = N_VLinearSumVectorArray(nvec, a, X_1d_ptr, b,
                                                   Y_1d_ptr, Z_1d_ptr);
      return lambda_result;
    };

    return N_VLinearSumVectorArray_adapt_arr_ptr_to_std_vector(nvec, a, X_1d, b,
                                                               Y_1d, Z_1d);
  },
  nb::arg("nvec"), nb::arg("a"), nb::arg("X_1d"), nb::arg("b"), nb::arg("Y_1d"),
  nb::arg("Z_1d"));

m.def(
  "N_VScaleVectorArray",
  [](int nvec, sundials4py::Array1d c_1d, std::vector<N_Vector> X_1d,
     std::vector<N_Vector> Z_1d) -> SUNErrCode
  {
    auto N_VScaleVectorArray_adapt_arr_ptr_to_std_vector =
      [](int nvec, sundials4py::Array1d c_1d, std::vector<N_Vector> X_1d,
         std::vector<N_Vector> Z_1d) -> SUNErrCode
    {
      sunrealtype* c_1d_ptr = reinterpret_cast<sunrealtype*>(c_1d.data());
      N_Vector* X_1d_ptr =
        reinterpret_cast<N_Vector*>(X_1d.empty() ? nullptr : X_1d.data());
      N_Vector* Z_1d_ptr =
        reinterpret_cast<N_Vector*>(Z_1d.empty() ? nullptr : Z_1d.data());

      auto lambda_result = N_VScaleVectorArray(nvec, c_1d_ptr, X_1d_ptr,
                                               Z_1d_ptr);
      return lambda_result;
    };

    return N_VScaleVectorArray_adapt_arr_ptr_to_std_vector(nvec, c_1d, X_1d,
                                                           Z_1d);
  },
  nb::arg("nvec"), nb::arg("c_1d"), nb::arg("X_1d"), nb::arg("Z_1d"));

m.def(
  "N_VConstVectorArray",
  [](int nvec, sunrealtype c, std::vector<N_Vector> Z_1d) -> SUNErrCode
  {
    auto N_VConstVectorArray_adapt_arr_ptr_to_std_vector =
      [](int nvec, sunrealtype c, std::vector<N_Vector> Z_1d) -> SUNErrCode
    {
      N_Vector* Z_1d_ptr =
        reinterpret_cast<N_Vector*>(Z_1d.empty() ? nullptr : Z_1d.data());

      auto lambda_result = N_VConstVectorArray(nvec, c, Z_1d_ptr);
      return lambda_result;
    };

    return N_VConstVectorArray_adapt_arr_ptr_to_std_vector(nvec, c, Z_1d);
  },
  nb::arg("nvec"), nb::arg("c"), nb::arg("Z_1d"));

m.def(
  "N_VWrmsNormVectorArray",
  [](int nvec, std::vector<N_Vector> X_1d, std::vector<N_Vector> W_1d,
     sundials4py::Array1d nrm_1d) -> SUNErrCode
  {
    auto N_VWrmsNormVectorArray_adapt_arr_ptr_to_std_vector =
      [](int nvec, std::vector<N_Vector> X_1d, std::vector<N_Vector> W_1d,
         sundials4py::Array1d nrm_1d) -> SUNErrCode
    {
      N_Vector* X_1d_ptr =
        reinterpret_cast<N_Vector*>(X_1d.empty() ? nullptr : X_1d.data());
      N_Vector* W_1d_ptr =
        reinterpret_cast<N_Vector*>(W_1d.empty() ? nullptr : W_1d.data());
      sunrealtype* nrm_1d_ptr = reinterpret_cast<sunrealtype*>(nrm_1d.data());

      auto lambda_result = N_VWrmsNormVectorArray(nvec, X_1d_ptr, W_1d_ptr,
                                                  nrm_1d_ptr);
      return lambda_result;
    };

    return N_VWrmsNormVectorArray_adapt_arr_ptr_to_std_vector(nvec, X_1d, W_1d,
                                                              nrm_1d);
  },
  nb::arg("nvec"), nb::arg("X_1d"), nb::arg("W_1d"), nb::arg("nrm_1d"));

m.def(
  "N_VWrmsNormMaskVectorArray",
  [](int nvec, std::vector<N_Vector> X_1d, std::vector<N_Vector> W_1d,
     N_Vector id, sundials4py::Array1d nrm_1d) -> SUNErrCode
  {
    auto N_VWrmsNormMaskVectorArray_adapt_arr_ptr_to_std_vector =
      [](int nvec, std::vector<N_Vector> X_1d, std::vector<N_Vector> W_1d,
         N_Vector id, sundials4py::Array1d nrm_1d) -> SUNErrCode
    {
      N_Vector* X_1d_ptr =
        reinterpret_cast<N_Vector*>(X_1d.empty() ? nullptr : X_1d.data());
      N_Vector* W_1d_ptr =
        reinterpret_cast<N_Vector*>(W_1d.empty() ? nullptr : W_1d.data());
      sunrealtype* nrm_1d_ptr = reinterpret_cast<sunrealtype*>(nrm_1d.data());

      auto lambda_result = N_VWrmsNormMaskVectorArray(nvec, X_1d_ptr, W_1d_ptr,
                                                      id, nrm_1d_ptr);
      return lambda_result;
    };

    return N_VWrmsNormMaskVectorArray_adapt_arr_ptr_to_std_vector(nvec, X_1d,
                                                                  W_1d, id,
                                                                  nrm_1d);
  },
  nb::arg("nvec"), nb::arg("X_1d"), nb::arg("W_1d"), nb::arg("id"),
  nb::arg("nrm_1d"));

m.def("N_VDotProdLocal", N_VDotProdLocal, nb::arg("x"), nb::arg("y"));

m.def("N_VMaxNormLocal", N_VMaxNormLocal, nb::arg("x"));

m.def("N_VMinLocal", N_VMinLocal, nb::arg("x"));

m.def("N_VL1NormLocal", N_VL1NormLocal, nb::arg("x"));

m.def("N_VWSqrSumLocal", N_VWSqrSumLocal, nb::arg("x"), nb::arg("w"));

m.def("N_VWSqrSumMaskLocal", N_VWSqrSumMaskLocal, nb::arg("x"), nb::arg("w"),
      nb::arg("id"));

m.def("N_VInvTestLocal", N_VInvTestLocal, nb::arg("x"), nb::arg("z"));

m.def("N_VConstrMaskLocal", N_VConstrMaskLocal, nb::arg("c"), nb::arg("x"),
      nb::arg("m"));

m.def("N_VMinQuotientLocal", N_VMinQuotientLocal, nb::arg("num"),
      nb::arg("denom"));

m.def(
  "N_VDotProdMultiLocal",
  [](int nvec, N_Vector x, std::vector<N_Vector> Y_1d,
     sundials4py::Array1d dotprods_1d) -> SUNErrCode
  {
    auto N_VDotProdMultiLocal_adapt_arr_ptr_to_std_vector =
      [](int nvec, N_Vector x, std::vector<N_Vector> Y_1d,
         sundials4py::Array1d dotprods_1d) -> SUNErrCode
    {
      N_Vector* Y_1d_ptr =
        reinterpret_cast<N_Vector*>(Y_1d.empty() ? nullptr : Y_1d.data());
      sunrealtype* dotprods_1d_ptr =
        reinterpret_cast<sunrealtype*>(dotprods_1d.data());

      auto lambda_result = N_VDotProdMultiLocal(nvec, x, Y_1d_ptr,
                                                dotprods_1d_ptr);
      return lambda_result;
    };

    return N_VDotProdMultiLocal_adapt_arr_ptr_to_std_vector(nvec, x, Y_1d,
                                                            dotprods_1d);
  },
  nb::arg("nvec"), nb::arg("x"), nb::arg("Y_1d"), nb::arg("dotprods_1d"));

m.def(
  "N_VDotProdMultiAllReduce",
  [](int nvec_total, N_Vector x, sundials4py::Array1d sum_1d) -> SUNErrCode
  {
    auto N_VDotProdMultiAllReduce_adapt_arr_ptr_to_std_vector =
      [](int nvec_total, N_Vector x, sundials4py::Array1d sum_1d) -> SUNErrCode
    {
      sunrealtype* sum_1d_ptr = reinterpret_cast<sunrealtype*>(sum_1d.data());

      auto lambda_result = N_VDotProdMultiAllReduce(nvec_total, x, sum_1d_ptr);
      return lambda_result;
    };

    return N_VDotProdMultiAllReduce_adapt_arr_ptr_to_std_vector(nvec_total, x,
                                                                sum_1d);
  },
  nb::arg("nvec_total"), nb::arg("x"), nb::arg("sum_1d"));

m.def("N_VPrint", N_VPrint, nb::arg("v"));

m.def("N_VPrintFile", N_VPrintFile, nb::arg("v"), nb::arg("outfile"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
