set(example cvAdvDiff_kry_hip)
sundials_add_executable(${example} ${example}.cpp)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
target_link_libraries(${example} PRIVATE hip::device sundials_cvode sundials_nvechip
                          ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(${example} PRIVATE sundials_cvode_fused_hip)
endif()
sundials_add_test(${example} ${example} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ${example}.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES ${example}.cpp ${example}.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/hip/cvAdvDiff_kry_hip)
endif()
