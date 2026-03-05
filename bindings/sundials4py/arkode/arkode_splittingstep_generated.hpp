// #ifndef ARKODE_SPLITTINGSTEP_H_
//
// #ifdef __cplusplus
// #endif
//

auto pyClassSplittingStepCoefficientsMem =
  nb::class_<SplittingStepCoefficientsMem>(m, "SplittingStepCoefficientsMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyEnumARKODE_SplittingCoefficientsID =
  nb::enum_<ARKODE_SplittingCoefficientsID>(m, "ARKODE_SplittingCoefficientsID",
                                            nb::is_arithmetic(), "")
    .value("ARKODE_SPLITTING_NONE", ARKODE_SPLITTING_NONE,
           "ensure enum is signed int")
    .value("ARKODE_SPLITTING_LIE_TROTTER_1_1_2",
           ARKODE_SPLITTING_LIE_TROTTER_1_1_2, "")
    .value("ARKODE_MIN_SPLITTING_NUM", ARKODE_MIN_SPLITTING_NUM, "")
    .value("ARKODE_SPLITTING_STRANG_2_2_2", ARKODE_SPLITTING_STRANG_2_2_2, "")
    .value("ARKODE_SPLITTING_BEST_2_2_2", ARKODE_SPLITTING_BEST_2_2_2, "")
    .value("ARKODE_SPLITTING_SUZUKI_3_3_2", ARKODE_SPLITTING_SUZUKI_3_3_2, "")
    .value("ARKODE_SPLITTING_RUTH_3_3_2", ARKODE_SPLITTING_RUTH_3_3_2, "")
    .value("ARKODE_SPLITTING_YOSHIDA_4_4_2", ARKODE_SPLITTING_YOSHIDA_4_4_2, "")
    .value("ARKODE_SPLITTING_YOSHIDA_8_6_2", ARKODE_SPLITTING_YOSHIDA_8_6_2, "")
    .value("ARKODE_MAX_SPLITTING_NUM", ARKODE_MAX_SPLITTING_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def(
  "SplittingStepCoefficients_Create",
  [](int sequential_methods, int stages, int partitions, int order,
     sundials4py::Array1d alpha_1d, sundials4py::Array1d beta_1d)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector =
      [](int sequential_methods, int stages, int partitions, int order,
         sundials4py::Array1d alpha_1d,
         sundials4py::Array1d beta_1d) -> SplittingStepCoefficients
    {
      sunrealtype* alpha_1d_ptr = reinterpret_cast<sunrealtype*>(alpha_1d.data());
      sunrealtype* beta_1d_ptr = reinterpret_cast<sunrealtype*>(beta_1d.data());

      auto lambda_result =
        SplittingStepCoefficients_Create(sequential_methods, stages, partitions,
                                         order, alpha_1d_ptr, beta_1d_ptr);
      return lambda_result;
    };
    auto SplittingStepCoefficients_Create_adapt_return_type_to_shared_ptr =
      [&SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector](int sequential_methods,
                                                                      int stages,
                                                                      int partitions,
                                                                      int order,
                                                                      sundials4py::Array1d
                                                                        alpha_1d,
                                                                      sundials4py::Array1d
                                                                        beta_1d)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result =
        SplittingStepCoefficients_Create_adapt_arr_ptr_to_std_vector(sequential_methods,
                                                                     stages,
                                                                     partitions,
                                                                     order,
                                                                     alpha_1d,
                                                                     beta_1d);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_Create_adapt_return_type_to_shared_ptr(sequential_methods,
                                                                            stages,
                                                                            partitions,
                                                                            order,
                                                                            alpha_1d,
                                                                            beta_1d);
  },
  nb::arg("sequential_methods"), nb::arg("stages"), nb::arg("partitions"),
  nb::arg("order"), nb::arg("alpha_1d"), nb::arg("beta_1d"));

m.def(
  "SplittingStepCoefficients_Copy",
  [](SplittingStepCoefficients coefficients)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_Copy_adapt_return_type_to_shared_ptr =
      [](SplittingStepCoefficients coefficients)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_Copy(coefficients);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_Copy_adapt_return_type_to_shared_ptr(
      coefficients);
  },
  nb::arg("coefficients"));

m.def("SplittingStepCoefficients_Write", SplittingStepCoefficients_Write,
      nb::arg("coefficients"), nb::arg("outfile"));

m.def(
  "SplittingStepCoefficients_LoadCoefficients",
  [](ARKODE_SplittingCoefficientsID id)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_LoadCoefficients_adapt_return_type_to_shared_ptr =
      [](ARKODE_SplittingCoefficientsID id)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_LoadCoefficients(id);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_LoadCoefficients_adapt_return_type_to_shared_ptr(
      id);
  },
  nb::arg("id"));

m.def(
  "SplittingStepCoefficients_LoadCoefficientsByName",
  [](const char* name)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_LoadCoefficientsByName_adapt_return_type_to_shared_ptr =
      [](const char* name)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_LoadCoefficientsByName(name);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_LoadCoefficientsByName_adapt_return_type_to_shared_ptr(
      name);
  },
  nb::arg("name"));

m.def("SplittingStepCoefficients_IDToName", SplittingStepCoefficients_IDToName,
      nb::arg("id"));

m.def(
  "SplittingStepCoefficients_LieTrotter",
  [](int partitions)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_LieTrotter_adapt_return_type_to_shared_ptr =
      [](int partitions)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_LieTrotter(partitions);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_LieTrotter_adapt_return_type_to_shared_ptr(
      partitions);
  },
  nb::arg("partitions"));

m.def(
  "SplittingStepCoefficients_Strang",
  [](int partitions)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_Strang_adapt_return_type_to_shared_ptr =
      [](int partitions)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_Strang(partitions);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_Strang_adapt_return_type_to_shared_ptr(
      partitions);
  },
  nb::arg("partitions"));

m.def(
  "SplittingStepCoefficients_Parallel",
  [](int partitions)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_Parallel_adapt_return_type_to_shared_ptr =
      [](int partitions)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_Parallel(partitions);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_Parallel_adapt_return_type_to_shared_ptr(
      partitions);
  },
  nb::arg("partitions"));

m.def(
  "SplittingStepCoefficients_SymmetricParallel",
  [](int partitions)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_SymmetricParallel_adapt_return_type_to_shared_ptr =
      [](int partitions)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_SymmetricParallel(partitions);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_SymmetricParallel_adapt_return_type_to_shared_ptr(
      partitions);
  },
  nb::arg("partitions"));

m.def(
  "SplittingStepCoefficients_ThirdOrderSuzuki",
  [](int partitions)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_ThirdOrderSuzuki_adapt_return_type_to_shared_ptr =
      [](int partitions)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_ThirdOrderSuzuki(partitions);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_ThirdOrderSuzuki_adapt_return_type_to_shared_ptr(
      partitions);
  },
  nb::arg("partitions"));

m.def(
  "SplittingStepCoefficients_TripleJump",
  [](int partitions, int order)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_TripleJump_adapt_return_type_to_shared_ptr =
      [](int partitions, int order)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_TripleJump(partitions,
                                                                order);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_TripleJump_adapt_return_type_to_shared_ptr(partitions,
                                                                                order);
  },
  nb::arg("partitions"), nb::arg("order"));

m.def(
  "SplittingStepCoefficients_SuzukiFractal",
  [](int partitions, int order)
    -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
  {
    auto SplittingStepCoefficients_SuzukiFractal_adapt_return_type_to_shared_ptr =
      [](int partitions, int order)
      -> std::shared_ptr<std::remove_pointer_t<SplittingStepCoefficients>>
    {
      auto lambda_result = SplittingStepCoefficients_SuzukiFractal(partitions,
                                                                   order);

      return our_make_shared<std::remove_pointer_t<SplittingStepCoefficients>,
                             SplittingStepCoefficientsDeleter>(lambda_result);
    };

    return SplittingStepCoefficients_SuzukiFractal_adapt_return_type_to_shared_ptr(partitions,
                                                                                   order);
  },
  nb::arg("partitions"), nb::arg("order"));

m.def("SplittingStepSetCoefficients", SplittingStepSetCoefficients,
      nb::arg("arkode_mem"), nb::arg("coefficients"));

m.def(
  "SplittingStepGetNumEvolves",
  [](void* arkode_mem, int partition) -> std::tuple<int, long>
  {
    auto SplittingStepGetNumEvolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int partition) -> std::tuple<int, long>
    {
      long evolves_adapt_modifiable;

      int r = SplittingStepGetNumEvolves(arkode_mem, partition,
                                         &evolves_adapt_modifiable);
      return std::make_tuple(r, evolves_adapt_modifiable);
    };

    return SplittingStepGetNumEvolves_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                           partition);
  },
  nb::arg("arkode_mem"), nb::arg("partition"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
