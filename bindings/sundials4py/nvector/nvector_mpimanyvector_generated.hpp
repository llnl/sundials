// #ifndef SUNDIALS_NVECTOR_MPIMANYVECTOR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_MPIManyVector =
  nb::class_<_N_VectorContent_MPIManyVector>(m, "_N_VectorContent_MPIManyVector",
                                             "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("N_VGetNumSubvectors_MPIManyVector", N_VGetNumSubvectors_MPIManyVector,
      nb::arg("v"));

m.def("N_VGetSubvectorLocalLength_MPIManyVector",
      N_VGetSubvectorLocalLength_MPIManyVector, nb::arg("v"), nb::arg("vec_num"));

m.def("N_VEnableFusedOps_MPIManyVector", N_VEnableFusedOps_MPIManyVector,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearCombination_MPIManyVector",
      N_VEnableLinearCombination_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMulti_MPIManyVector",
      N_VEnableScaleAddMulti_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableDotProdMulti_MPIManyVector",
      N_VEnableDotProdMulti_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearSumVectorArray_MPIManyVector",
      N_VEnableLinearSumVectorArray_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleVectorArray_MPIManyVector",
      N_VEnableScaleVectorArray_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableConstVectorArray_MPIManyVector",
      N_VEnableConstVectorArray_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormVectorArray_MPIManyVector",
      N_VEnableWrmsNormVectorArray_MPIManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormMaskVectorArray_MPIManyVector",
      N_VEnableWrmsNormMaskVectorArray_MPIManyVector, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableDotProdMultiLocal_MPIManyVector",
      N_VEnableDotProdMultiLocal_MPIManyVector, nb::arg("v"), nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
