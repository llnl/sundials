# CVODE RAJA C++ example
set(CVODE_examples "cvAdvDiff_kry_raja.cpp\;develop")
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  set(_fused_stubs sundials_cvode_fused_stubs)
endif()
if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(_other_libs OpenMP::OpenMP_CXX)
endif()
set_source_files_properties(cvAdvDiff_kry_raja.cpp PROPERTIES LANGUAGE ${_lang})
sundials_add_executable(cvAdvDiff_kry_raja cvAdvDiff_kry_raja.cpp)
set_target_properties(cvAdvDiff_kry_raja PROPERTIES FOLDER "Examples")
target_include_directories(cvAdvDiff_kry_raja PRIVATE .)
target_link_libraries(cvAdvDiff_kry_raja
  PRIVATE sundials_cvode sundials_nvecraja ${_fused_stubs} ${EXE_EXTRA_LINK_LIBS}
          ${_other_libs})
sundials_add_test(cvAdvDiff_kry_raja cvAdvDiff_kry_raja
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_kry_raja.out EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_kry_raja.cpp cvAdvDiff_kry_raja.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cvAdvDiff_kry_raja)
endif()
