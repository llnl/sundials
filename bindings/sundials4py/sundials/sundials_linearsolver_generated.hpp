// #ifndef _SUNDIALS_ITERATIVE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNPrecType = nb::enum_<SUNPrecType>(m, "SUNPrecType",
                                                nb::is_arithmetic(), "")
                           .value("SUN_PREC_NONE", SUN_PREC_NONE, "")
                           .value("SUN_PREC_LEFT", SUN_PREC_LEFT, "")
                           .value("SUN_PREC_RIGHT", SUN_PREC_RIGHT, "")
                           .value("SUN_PREC_BOTH", SUN_PREC_BOTH, "")
                           .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyEnumSUNGramSchmidtType =
  nb::enum_<SUNGramSchmidtType>(m, "SUNGramSchmidtType", nb::is_arithmetic(), "")
    .value("SUN_MODIFIED_GS", SUN_MODIFIED_GS, "")
    .value("SUN_CLASSICAL_GS", SUN_CLASSICAL_GS, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def(
  "SUNModifiedGS",
  [](std::vector<N_Vector> v_1d, sundials4py::Array1d h_2d, int k,
     int p) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNModifiedGS_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> v_1d, sundials4py::Array1d h_2d, int k, int p,
         sunrealtype* new_vk_norm) -> SUNErrCode
    {
      N_Vector* v_1d_ptr     = v_1d.empty() ? nullptr : v_1d.data();
      sunrealtype** h_2d_ptr = reinterpret_cast<sunrealtype**>(
        h_2d.size() == 0 ? nullptr : h_2d.data());

      auto lambda_result = SUNModifiedGS(v_1d_ptr, h_2d_ptr, k, p, new_vk_norm);
      return lambda_result;
    };
    auto SUNModifiedGS_adapt_modifiable_immutable_to_return =
      [&SUNModifiedGS_adapt_arr_ptr_to_std_vector](std::vector<N_Vector> v_1d,
                                                   sundials4py::Array1d h_2d,
                                                   int k, int p)
      -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype new_vk_norm_adapt_modifiable;

      SUNErrCode r =
        SUNModifiedGS_adapt_arr_ptr_to_std_vector(v_1d, h_2d, k, p,
                                                  &new_vk_norm_adapt_modifiable);
      return std::make_tuple(r, new_vk_norm_adapt_modifiable);
    };

    return SUNModifiedGS_adapt_modifiable_immutable_to_return(v_1d, h_2d, k, p);
  },
  nb::arg("v_1d"), nb::arg("h_2d"), nb::arg("k"), nb::arg("p"));

m.def(
  "SUNClassicalGS",
  [](std::vector<N_Vector> v_1d, sundials4py::Array1d h_2d, int k, int p,
     sundials4py::Array1d stemp_1d,
     std::vector<N_Vector> vtemp_1d) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNClassicalGS_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> v_1d, sundials4py::Array1d h_2d, int k, int p,
         sunrealtype* new_vk_norm, sundials4py::Array1d stemp_1d,
         std::vector<N_Vector> vtemp_1d) -> SUNErrCode
    {
      N_Vector* v_1d_ptr     = v_1d.empty() ? nullptr : v_1d.data();
      sunrealtype** h_2d_ptr = reinterpret_cast<sunrealtype**>(
        h_2d.size() == 0 ? nullptr : h_2d.data());
      sunrealtype* stemp_1d_ptr = stemp_1d.size() == 0 ? nullptr
                                                       : stemp_1d.data();
      N_Vector* vtemp_1d_ptr    = vtemp_1d.empty() ? nullptr : vtemp_1d.data();

      auto lambda_result = SUNClassicalGS(v_1d_ptr, h_2d_ptr, k, p, new_vk_norm,
                                          stemp_1d_ptr, vtemp_1d_ptr);
      return lambda_result;
    };
    auto SUNClassicalGS_adapt_modifiable_immutable_to_return =
      [&SUNClassicalGS_adapt_arr_ptr_to_std_vector](std::vector<N_Vector> v_1d,
                                                    sundials4py::Array1d h_2d,
                                                    int k, int p,
                                                    sundials4py::Array1d stemp_1d,
                                                    std::vector<N_Vector> vtemp_1d)
      -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype new_vk_norm_adapt_modifiable;

      SUNErrCode r =
        SUNClassicalGS_adapt_arr_ptr_to_std_vector(v_1d, h_2d, k, p,
                                                   &new_vk_norm_adapt_modifiable,
                                                   stemp_1d, vtemp_1d);
      return std::make_tuple(r, new_vk_norm_adapt_modifiable);
    };

    return SUNClassicalGS_adapt_modifiable_immutable_to_return(v_1d, h_2d, k, p,
                                                               stemp_1d,
                                                               vtemp_1d);
  },
  nb::arg("v_1d"), nb::arg("h_2d"), nb::arg("k"), nb::arg("p"),
  nb::arg("stemp_1d"), nb::arg("vtemp_1d"));

m.def(
  "SUNQRfact",
  [](int n, sundials4py::Array1d h_2d, sundials4py::Array1d q_1d, int job) -> int
  {
    auto SUNQRfact_adapt_arr_ptr_to_std_vector =
      [](int n, sundials4py::Array1d h_2d, sundials4py::Array1d q_1d,
         int job) -> int
    {
      sunrealtype** h_2d_ptr = reinterpret_cast<sunrealtype**>(
        h_2d.size() == 0 ? nullptr : h_2d.data());
      sunrealtype* q_1d_ptr = q_1d.size() == 0 ? nullptr : q_1d.data();

      auto lambda_result = SUNQRfact(n, h_2d_ptr, q_1d_ptr, job);
      return lambda_result;
    };

    return SUNQRfact_adapt_arr_ptr_to_std_vector(n, h_2d, q_1d, job);
  },
  nb::arg("n"), nb::arg("h_2d"), nb::arg("q_1d"), nb::arg("job"));

m.def(
  "SUNQRsol",
  [](int n, sundials4py::Array1d h_2d, sundials4py::Array1d q_1d,
     sundials4py::Array1d b_1d) -> int
  {
    auto SUNQRsol_adapt_arr_ptr_to_std_vector =
      [](int n, sundials4py::Array1d h_2d, sundials4py::Array1d q_1d,
         sundials4py::Array1d b_1d) -> int
    {
      sunrealtype** h_2d_ptr = reinterpret_cast<sunrealtype**>(
        h_2d.size() == 0 ? nullptr : h_2d.data());
      sunrealtype* q_1d_ptr = q_1d.size() == 0 ? nullptr : q_1d.data();
      sunrealtype* b_1d_ptr = b_1d.size() == 0 ? nullptr : b_1d.data();

      auto lambda_result = SUNQRsol(n, h_2d_ptr, q_1d_ptr, b_1d_ptr);
      return lambda_result;
    };

    return SUNQRsol_adapt_arr_ptr_to_std_vector(n, h_2d, q_1d, b_1d);
  },
  nb::arg("n"), nb::arg("h_2d"), nb::arg("q_1d"), nb::arg("b_1d"));

m.def(
  "SUNQRAdd_MGS",
  [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_MGS_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr    = Q_1d.empty() ? nullptr : Q_1d.data();
      sunrealtype* R_1d_ptr = R_1d.size() == 0 ? nullptr : R_1d.data();

      auto lambda_result = SUNQRAdd_MGS(Q_1d_ptr, R_1d_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_MGS_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                    QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_ICWY",
  [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_ICWY_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr    = Q_1d.empty() ? nullptr : Q_1d.data();
      sunrealtype* R_1d_ptr = R_1d.size() == 0 ? nullptr : R_1d.data();

      auto lambda_result = SUNQRAdd_ICWY(Q_1d_ptr, R_1d_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_ICWY_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                     QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_ICWY_SB",
  [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_ICWY_SB_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr    = Q_1d.empty() ? nullptr : Q_1d.data();
      sunrealtype* R_1d_ptr = R_1d.size() == 0 ? nullptr : R_1d.data();

      auto lambda_result = SUNQRAdd_ICWY_SB(Q_1d_ptr, R_1d_ptr, df, m, mMax,
                                            QRdata);
      return lambda_result;
    };

    return SUNQRAdd_ICWY_SB_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                        QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_CGS2",
  [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_CGS2_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr    = Q_1d.empty() ? nullptr : Q_1d.data();
      sunrealtype* R_1d_ptr = R_1d.size() == 0 ? nullptr : R_1d.data();

      auto lambda_result = SUNQRAdd_CGS2(Q_1d_ptr, R_1d_ptr, df, m, mMax, QRdata);
      return lambda_result;
    };

    return SUNQRAdd_CGS2_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                     QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_DCGS2",
  [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_DCGS2_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr    = Q_1d.empty() ? nullptr : Q_1d.data();
      sunrealtype* R_1d_ptr = R_1d.size() == 0 ? nullptr : R_1d.data();

      auto lambda_result = SUNQRAdd_DCGS2(Q_1d_ptr, R_1d_ptr, df, m, mMax,
                                          QRdata);
      return lambda_result;
    };

    return SUNQRAdd_DCGS2_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m, mMax,
                                                      QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));

m.def(
  "SUNQRAdd_DCGS2_SB",
  [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df, int m,
     int mMax, void* QRdata) -> SUNErrCode
  {
    auto SUNQRAdd_DCGS2_SB_adapt_arr_ptr_to_std_vector =
      [](std::vector<N_Vector> Q_1d, sundials4py::Array1d R_1d, N_Vector df,
         int m, int mMax, void* QRdata) -> SUNErrCode
    {
      N_Vector* Q_1d_ptr    = Q_1d.empty() ? nullptr : Q_1d.data();
      sunrealtype* R_1d_ptr = R_1d.size() == 0 ? nullptr : R_1d.data();

      auto lambda_result = SUNQRAdd_DCGS2_SB(Q_1d_ptr, R_1d_ptr, df, m, mMax,
                                             QRdata);
      return lambda_result;
    };

    return SUNQRAdd_DCGS2_SB_adapt_arr_ptr_to_std_vector(Q_1d, R_1d, df, m,
                                                         mMax, QRdata);
  },
  nb::arg("Q_1d"), nb::arg("R_1d"), nb::arg("df"), nb::arg("m"),
  nb::arg("mMax"), nb::arg("QRdata"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _SUNLINEARSOLVER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNLinearSolver_Type =
  nb::enum_<SUNLinearSolver_Type>(m, "SUNLinearSolver_Type",
                                  nb::is_arithmetic(), "")
    .value("SUNLINEARSOLVER_DIRECT", SUNLINEARSOLVER_DIRECT, "")
    .value("SUNLINEARSOLVER_ITERATIVE", SUNLINEARSOLVER_ITERATIVE, "")
    .value("SUNLINEARSOLVER_MATRIX_ITERATIVE", SUNLINEARSOLVER_MATRIX_ITERATIVE,
           "")
    .value("SUNLINEARSOLVER_MATRIX_EMBEDDED", SUNLINEARSOLVER_MATRIX_EMBEDDED, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyEnumSUNLinearSolver_ID =
  nb::enum_<SUNLinearSolver_ID>(m, "SUNLinearSolver_ID", nb::is_arithmetic(), "")
    .value("SUNLINEARSOLVER_BAND", SUNLINEARSOLVER_BAND, "")
    .value("SUNLINEARSOLVER_DENSE", SUNLINEARSOLVER_DENSE, "")
    .value("SUNLINEARSOLVER_KLU", SUNLINEARSOLVER_KLU, "")
    .value("SUNLINEARSOLVER_LAPACKBAND", SUNLINEARSOLVER_LAPACKBAND, "")
    .value("SUNLINEARSOLVER_LAPACKDENSE", SUNLINEARSOLVER_LAPACKDENSE, "")
    .value("SUNLINEARSOLVER_PCG", SUNLINEARSOLVER_PCG, "")
    .value("SUNLINEARSOLVER_SPBCGS", SUNLINEARSOLVER_SPBCGS, "")
    .value("SUNLINEARSOLVER_SPFGMR", SUNLINEARSOLVER_SPFGMR, "")
    .value("SUNLINEARSOLVER_SPGMR", SUNLINEARSOLVER_SPGMR, "")
    .value("SUNLINEARSOLVER_SPTFQMR", SUNLINEARSOLVER_SPTFQMR, "")
    .value("SUNLINEARSOLVER_SUPERLUDIST", SUNLINEARSOLVER_SUPERLUDIST, "")
    .value("SUNLINEARSOLVER_SUPERLUMT", SUNLINEARSOLVER_SUPERLUMT, "")
    .value("SUNLINEARSOLVER_CUSOLVERSP_BATCHQR",
           SUNLINEARSOLVER_CUSOLVERSP_BATCHQR, "")
    .value("SUNLINEARSOLVER_MAGMADENSE", SUNLINEARSOLVER_MAGMADENSE, "")
    .value("SUNLINEARSOLVER_ONEMKLDENSE", SUNLINEARSOLVER_ONEMKLDENSE, "")
    .value("SUNLINEARSOLVER_GINKGO", SUNLINEARSOLVER_GINKGO, "")
    .value("SUNLINEARSOLVER_GINKGOBATCH", SUNLINEARSOLVER_GINKGOBATCH, "")
    .value("SUNLINEARSOLVER_KOKKOSDENSE", SUNLINEARSOLVER_KOKKOSDENSE, "")
    .value("SUNLINEARSOLVER_CUSTOM", SUNLINEARSOLVER_CUSTOM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNLinearSolver_Ops =
  nb::class_<_generic_SUNLinearSolver_Ops>(m,
                                           "_generic_SUNLinearSolver_Ops", "Structure containing function pointers to linear solver operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNLinearSolver =
  nb::class_<_generic_SUNLinearSolver>(m,
                                       "_generic_SUNLinearSolver", " A linear solver is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of linear solver\n   operations corresponding to that implementation.")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNLinSolGetType", SUNLinSolGetType, nb::arg("S"));

m.def("SUNLinSolGetID", SUNLinSolGetID, nb::arg("S"));

m.def("SUNLinSolSetScalingVectors", SUNLinSolSetScalingVectors, nb::arg("S"),
      nb::arg("s1"), nb::arg("s2"));

m.def("SUNLinSolSetZeroGuess", SUNLinSolSetZeroGuess, nb::arg("S"),
      nb::arg("onoff"));

m.def("SUNLinSolInitialize", SUNLinSolInitialize, nb::arg("S"));

m.def(
  "SUNLinSolSetup",
  [](SUNLinearSolver S, std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto SUNLinSolSetup_adapt_optional_arg_with_default_null =
      [](SUNLinearSolver S, std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = SUNLinSolSetup(S, A_adapt_default_null);
      return lambda_result;
    };

    return SUNLinSolSetup_adapt_optional_arg_with_default_null(S, A);
  },
  nb::arg("S"), nb::arg("A").none() = nb::none());

m.def("SUNLinSolNumIters", SUNLinSolNumIters, nb::arg("S"));

m.def("SUNLinSolResNorm", SUNLinSolResNorm, nb::arg("S"));

m.def("SUNLinSolResid", SUNLinSolResid, nb::arg("S"),
      "nb::rv_policy::reference", nb::rv_policy::reference);

m.def("SUNLinSolLastFlag", SUNLinSolLastFlag, nb::arg("S"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
