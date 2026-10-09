// #ifndef _NVECTOR_MANY_VECTOR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_ManyVector =
  nb::class_<_N_VectorContent_ManyVector>(m, "_N_VectorContent_ManyVector", "")
    .def(nb::init<>()) // implicit default constructor
  ;

sundials4py::scoped_def(m, "N_VGetSubvector_ManyVector",
                        N_VGetSubvector_ManyVector, nb::arg("v"),
                        nb::arg("vec_num"), "nb::rv_policy::reference",
                        nb::rv_policy::reference);

sundials4py::scoped_def(m, "N_VGetNumSubvectors_ManyVector",
                        N_VGetNumSubvectors_ManyVector, nb::arg("v"));

sundials4py::scoped_def(m, "N_VGetSubvectorLocalLength_ManyVector",
                        N_VGetSubvectorLocalLength_ManyVector, nb::arg("v"),
                        nb::arg("vec_num"));

sundials4py::scoped_def(m, "N_VEnableFusedOps_ManyVector",
                        N_VEnableFusedOps_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableLinearCombination_ManyVector",
                        N_VEnableLinearCombination_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableScaleAddMulti_ManyVector",
                        N_VEnableScaleAddMulti_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableDotProdMulti_ManyVector",
                        N_VEnableDotProdMulti_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableLinearSumVectorArray_ManyVector",
                        N_VEnableLinearSumVectorArray_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableScaleVectorArray_ManyVector",
                        N_VEnableScaleVectorArray_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableConstVectorArray_ManyVector",
                        N_VEnableConstVectorArray_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableWrmsNormVectorArray_ManyVector",
                        N_VEnableWrmsNormVectorArray_ManyVector, nb::arg("v"),
                        nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableWrmsNormMaskVectorArray_ManyVector",
                        N_VEnableWrmsNormMaskVectorArray_ManyVector,
                        nb::arg("v"), nb::arg("tf"));

sundials4py::scoped_def(m, "N_VEnableDotProdMultiLocal_ManyVector",
                        N_VEnableDotProdMultiLocal_ManyVector, nb::arg("v"),
                        nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
