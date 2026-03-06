// #ifndef _SUNDIALS_MEMORY_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNMemoryType = nb::enum_<SUNMemoryType>(m, "SUNMemoryType",
                                                    nb::is_arithmetic(), "")
                             .value("SUNMEMTYPE_HOST", SUNMEMTYPE_HOST,
                                    "pageable memory accessible on the host")
                             .value("SUNMEMTYPE_PINNED", SUNMEMTYPE_PINNED,
                                    "page-locked memory accessible on the host")
                             .value("SUNMEMTYPE_DEVICE", SUNMEMTYPE_DEVICE,
                                    "memory accessible from the device")
                             .value("SUNMEMTYPE_UVM", SUNMEMTYPE_UVM,
                                    "memory accessible from the host or device")
                             .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClassSUNMemoryHelper_ =
  nb::class_<SUNMemoryHelper_>(m, "SUNMemoryHelper_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClassSUNMemoryHelper_Ops_ =
  nb::class_<SUNMemoryHelper_Ops_>(m, "SUNMemoryHelper_Ops_", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "SUNMemoryHelper_Clone",
  [](SUNMemoryHelper param_0)
    -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
  {
    auto SUNMemoryHelper_Clone_adapt_return_type_to_shared_ptr =
      [](SUNMemoryHelper param_0)
      -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
    {
      auto lambda_result = SUNMemoryHelper_Clone(param_0);

      return our_make_shared<std::remove_pointer_t<SUNMemoryHelper>,
                             SUNMemoryHelperDeleter>(lambda_result);
    };

    return SUNMemoryHelper_Clone_adapt_return_type_to_shared_ptr(param_0);
  },
  nb::arg("param_0"));

m.def("SUNMemoryHelper_SetDefaultQueue", SUNMemoryHelper_SetDefaultQueue,
      nb::arg("param_0"), nb::arg("queue"));

m.def("SUNMemoryHelper_ImplementsRequiredOps",
      SUNMemoryHelper_ImplementsRequiredOps, nb::arg("param_0"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
