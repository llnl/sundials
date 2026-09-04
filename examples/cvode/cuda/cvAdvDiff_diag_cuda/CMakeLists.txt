set(example cvAdvDiff_diag_cuda)

set_source_files_properties(${example}.cu PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
sundials_add_executable(${example} ${example}.cu)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
target_link_libraries(${example} PRIVATE sundials_cvode sundials_nveccuda
                          ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(${example} PRIVATE sundials_cvode_fused_cuda)
endif()

foreach(example_args "0 0" "0 1" "1 1")
  string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  sundials_add_test(${test_name} ${example} TEST_ARGS ${example_args}
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ${test_name}.out EXAMPLE_TYPE develop)
endforeach()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out ${example}*.out)
  install(FILES ${example}.cu ${example_out}
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cuda/cvAdvDiff_diag_cuda)
endif()
