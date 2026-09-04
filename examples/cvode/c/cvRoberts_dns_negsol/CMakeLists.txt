sundials_add_executable(cvRoberts_dns_negsol cvRoberts_dns_negsol.c)
set_target_properties(cvRoberts_dns_negsol PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_dns_negsol PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_dns_negsol cvRoberts_dns_negsol
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_dns_negsol.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dns_negsol.c cvRoberts_dns_negsol.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_dns_negsol)
endif()
