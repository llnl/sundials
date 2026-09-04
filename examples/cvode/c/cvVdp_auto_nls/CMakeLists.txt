sundials_add_executable(cvVdp_auto_nls cvVdp_auto_nls.c)
set_target_properties(cvVdp_auto_nls PROPERTIES FOLDER "Examples")
target_link_libraries(cvVdp_auto_nls PRIVATE ${CVODE_C_LIBS}
                      sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})

sundials_add_test(cvVdp_auto_nls cvVdp_auto_nls
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvVdp_auto_nls.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvVdp_auto_nls.c cvVdp_auto_nls.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvVdp_auto_nls)
endif()
