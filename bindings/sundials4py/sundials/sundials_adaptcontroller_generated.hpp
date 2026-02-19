// #ifndef _SUNDIALS_ADAPTCONTROLLER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNAdaptController_Type =
  nb::enum_<SUNAdaptController_Type>(m, "SUNAdaptController_Type",
                                     nb::is_arithmetic(), "")
    .value("SUN_ADAPTCONTROLLER_NONE", SUN_ADAPTCONTROLLER_NONE, "")
    .value("SUN_ADAPTCONTROLLER_H", SUN_ADAPTCONTROLLER_H, "")
    .value("SUN_ADAPTCONTROLLER_MRI_H_TOL", SUN_ADAPTCONTROLLER_MRI_H_TOL, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNAdaptController_Ops =
  nb::class_<_generic_SUNAdaptController_Ops>(m,
                                              "_generic_SUNAdaptController_Ops", "Structure containing function pointers to controller operations")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNAdaptController =
  nb::class_<_generic_SUNAdaptController>(m,
                                          "_generic_SUNAdaptController", " A SUNAdaptController is a structure with an implementation-dependent\n   'content' field, and a pointer to a structure of\n   operations corresponding to that implementation.")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNAdaptController_GetType", SUNAdaptController_GetType, nb::arg("C"));

m.def(
  "SUNAdaptController_EstimateStep",
  [](SUNAdaptController C, sunrealtype h, int p,
     sunrealtype dsm) -> std::tuple<SUNErrCode, sunrealtype>
  {
    auto SUNAdaptController_EstimateStep_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C, sunrealtype h, int p,
         sunrealtype dsm) -> std::tuple<SUNErrCode, sunrealtype>
    {
      sunrealtype hnew_adapt_modifiable;

      SUNErrCode r = SUNAdaptController_EstimateStep(C, h, p, dsm,
                                                     &hnew_adapt_modifiable);
      return std::make_tuple(r, hnew_adapt_modifiable);
    };

    return SUNAdaptController_EstimateStep_adapt_modifiable_immutable_to_return(C,
                                                                                h,
                                                                                p,
                                                                                dsm);
  },
  nb::arg("C"), nb::arg("h"), nb::arg("p"), nb::arg("dsm"));

m.def(
  "SUNAdaptController_EstimateStepTol",
  [](SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
     sunrealtype DSM,
     sunrealtype dsm) -> std::tuple<SUNErrCode, sunrealtype, sunrealtype>
  {
    auto SUNAdaptController_EstimateStepTol_adapt_modifiable_immutable_to_return =
      [](SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
         sunrealtype DSM,
         sunrealtype dsm) -> std::tuple<SUNErrCode, sunrealtype, sunrealtype>
    {
      sunrealtype Hnew_adapt_modifiable;
      sunrealtype tolfacnew_adapt_modifiable;

      SUNErrCode r =
        SUNAdaptController_EstimateStepTol(C, H, tolfac, P, DSM, dsm,
                                           &Hnew_adapt_modifiable,
                                           &tolfacnew_adapt_modifiable);
      return std::make_tuple(r, Hnew_adapt_modifiable,
                             tolfacnew_adapt_modifiable);
    };

    return SUNAdaptController_EstimateStepTol_adapt_modifiable_immutable_to_return(C,
                                                                                   H,
                                                                                   tolfac,
                                                                                   P,
                                                                                   DSM,
                                                                                   dsm);
  },
  nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("P"), nb::arg("DSM"),
  nb::arg("dsm"));

m.def("SUNAdaptController_Reset", SUNAdaptController_Reset, nb::arg("C"));

m.def("SUNAdaptController_SetDefaults", SUNAdaptController_SetDefaults,
      nb::arg("C"));

m.def("SUNAdaptController_Write", SUNAdaptController_Write, nb::arg("C"),
      nb::arg("fptr"));

m.def("SUNAdaptController_SetErrorBias", SUNAdaptController_SetErrorBias,
      nb::arg("C"), nb::arg("bias"));

m.def("SUNAdaptController_UpdateH", SUNAdaptController_UpdateH, nb::arg("C"),
      nb::arg("h"), nb::arg("dsm"));

m.def("SUNAdaptController_UpdateMRIHTol", SUNAdaptController_UpdateMRIHTol,
      nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("DSM"),
      nb::arg("dsm"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
