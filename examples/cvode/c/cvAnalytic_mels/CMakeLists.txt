sundials_add_executable(cvAnalytic_mels cvAnalytic_mels.c)
set_target_properties(cvAnalytic_mels PROPERTIES FOLDER "Examples")
target_link_libraries(cvAnalytic_mels PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAnalytic_mels cvAnalytic_mels ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAnalytic_mels.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAnalytic_mels.c cvAnalytic_mels.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAnalytic_mels)
endif()
