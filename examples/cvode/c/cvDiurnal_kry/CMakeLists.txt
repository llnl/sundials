sundials_add_executable(cvDiurnal_kry cvDiurnal_kry.c)
set_target_properties(cvDiurnal_kry PROPERTIES FOLDER "Examples")
target_link_libraries(cvDiurnal_kry PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDiurnal_kry cvDiurnal_kry ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDiurnal_kry.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDiurnal_kry.c cvDiurnal_kry.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvDiurnal_kry)
endif()
