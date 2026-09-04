sundials_add_executable(cvRoberts_dns_constraints cvRoberts_dns_constraints.c)
set_target_properties(cvRoberts_dns_constraints PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_dns_constraints PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_dns_constraints cvRoberts_dns_constraints
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_dns_constraints.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dns_constraints.c cvRoberts_dns_constraints.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_dns_constraints)
endif()
