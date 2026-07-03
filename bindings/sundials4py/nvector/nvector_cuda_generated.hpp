// #ifndef _NVECTOR_CUDA_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_Cuda =
  nb::class_<_N_VectorContent_Cuda>(m, "_N_VectorContent_Cuda", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("N_VIsManagedMemory_Cuda", N_VIsManagedMemory_Cuda, nb::arg("x"));

m.def("N_VCopyToDevice_Cuda", N_VCopyToDevice_Cuda, nb::arg("v"));

m.def("N_VCopyFromDevice_Cuda", N_VCopyFromDevice_Cuda, nb::arg("v"));

m.def("N_VEnableFusedOps_Cuda", N_VEnableFusedOps_Cuda, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearCombination_Cuda", N_VEnableLinearCombination_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMulti_Cuda", N_VEnableScaleAddMulti_Cuda, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableDotProdMulti_Cuda", N_VEnableDotProdMulti_Cuda, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearSumVectorArray_Cuda", N_VEnableLinearSumVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleVectorArray_Cuda", N_VEnableScaleVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableConstVectorArray_Cuda", N_VEnableConstVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormVectorArray_Cuda", N_VEnableWrmsNormVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormMaskVectorArray_Cuda",
      N_VEnableWrmsNormMaskVectorArray_Cuda, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMultiVectorArray_Cuda",
      N_VEnableScaleAddMultiVectorArray_Cuda, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearCombinationVectorArray_Cuda",
      N_VEnableLinearCombinationVectorArray_Cuda, nb::arg("v"), nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
