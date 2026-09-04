# CVODE SuperLU_DIST C++ example
set(CVODE_examples "cvAdvDiff_sludist.cpp\;2 1\;2\;develop")
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  set(_fused_lib sundials_cvode_fused_stubs)
endif()
sundials_add_executable(cvAdvDiff_sludist cvAdvDiff_sludist.cpp)
set_target_properties(cvAdvDiff_sludist PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_sludist
  PRIVATE MPI::MPI_CXX sundials_cvode sundials_nvecparallel
          sundials_sunlinsolsuperludist ${_fused_lib} ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_INDEX_SIZE STREQUAL "32")
  set(_integer_precision 15)
  set(_float_precision 3)
else()
  set(_integer_precision default)
  set(_float_precision default)
endif()
sundials_add_test(cvAdvDiff_sludist cvAdvDiff_sludist TEST_ARGS 2 1 MPI_NPROCS 2
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_sludist_${SUNDIALS_INDEX_SIZE}.out
                  EXAMPLE_TYPE develop INTEGER_PRECISION ${_integer_precision}
                  FLOAT_PRECISION ${_float_precision})

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_sludist.cpp cvAdvDiff_sludist_32.out
                cvAdvDiff_sludist_64.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cvAdvDiff_sludist)
endif()
