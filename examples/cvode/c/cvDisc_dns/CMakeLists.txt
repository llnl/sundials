sundials_add_executable(cvDisc_dns cvDisc_dns.c)
set_target_properties(cvDisc_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cvDisc_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDisc_dns cvDisc_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDisc_dns.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDisc_dns.c cvDisc_dns.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvDisc_dns)
endif()
