// #ifndef SUNDIALS_NVECTOR_HIP_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_Hip =
  nb::class_<_N_VectorContent_Hip>(m, "_N_VectorContent_Hip", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("N_VIsManagedMemory_Hip", N_VIsManagedMemory_Hip, nb::arg("x"));

m.def("N_VEnableFusedOps_Hip", N_VEnableFusedOps_Hip, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearCombination_Hip", N_VEnableLinearCombination_Hip,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMulti_Hip", N_VEnableScaleAddMulti_Hip, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableDotProdMulti_Hip", N_VEnableDotProdMulti_Hip, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearSumVectorArray_Hip", N_VEnableLinearSumVectorArray_Hip,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleVectorArray_Hip", N_VEnableScaleVectorArray_Hip,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableConstVectorArray_Hip", N_VEnableConstVectorArray_Hip,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormVectorArray_Hip", N_VEnableWrmsNormVectorArray_Hip,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormMaskVectorArray_Hip",
      N_VEnableWrmsNormMaskVectorArray_Hip, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMultiVectorArray_Hip",
      N_VEnableScaleAddMultiVectorArray_Hip, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearCombinationVectorArray_Hip",
      N_VEnableLinearCombinationVectorArray_Hip, nb::arg("v"), nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
