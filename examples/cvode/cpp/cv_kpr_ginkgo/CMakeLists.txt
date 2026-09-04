# CVODE Ginkgo KPR example
set(cpu_examples "cv_kpr_ginkgo.cpp\;\;develop")
sundials_add_examples_ginkgo(cpu_examples TARGETS sundials_cvode BACKENDS REF OMP)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out cv_kpr_ginkgo*.out)
  install(FILES cv_kpr_ginkgo.cpp cv_kpr_ginkgo.hpp ${example_out}
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cv_kpr_ginkgo)
endif()
