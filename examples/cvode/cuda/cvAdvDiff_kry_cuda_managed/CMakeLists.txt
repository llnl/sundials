set(example cvAdvDiff_kry_cuda_managed)

set_source_files_properties(${example}.cu PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
sundials_add_executable(${example} ${example}.cu)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
target_link_libraries(${example} PRIVATE sundials_cvode sundials_nveccuda
                          ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(${example} PRIVATE sundials_cvode_fused_cuda)
endif()

sundials_add_test(${example} ${example} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ${example}.out EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES ${example}.cu ${example}.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cuda/cvAdvDiff_kry_cuda_managed)
endif()
