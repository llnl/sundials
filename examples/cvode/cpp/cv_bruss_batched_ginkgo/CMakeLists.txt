# CVODE Ginkgo batched Brusselator example
set(cpu_gpu_examples "cv_bruss_batched_ginkgo.cpp\;\;develop")
sundials_add_examples_ginkgo(cpu_gpu_examples TARGETS sundials_cvode
                              BACKENDS REF OMP CUDA HIP SYCL)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out cv_bruss_batched_ginkgo*.out)
  install(FILES cv_bruss_batched_ginkgo.cpp cv_bruss_batched_ginkgo.hpp ${example_out}
          "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cv_bruss_batched_ginkgo)
endif()
