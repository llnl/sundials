sundials_add_executable(cvRocket_dns cvRocket_dns.c)
set_target_properties(cvRocket_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cvRocket_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRocket_dns cvRocket_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRocket_dns.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRocket_dns.c cvRocket_dns.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRocket_dns)
endif()
