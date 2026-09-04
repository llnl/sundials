sundials_add_executable(cvPendulum_dns cvPendulum_dns.c)
set_target_properties(cvPendulum_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cvPendulum_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvPendulum_dns cvPendulum_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvPendulum_dns.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvPendulum_dns.c cvPendulum_dns.out plot_cvPendulum.py
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvPendulum_dns)
endif()
