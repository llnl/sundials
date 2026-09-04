sundials_add_executable(cvRoberts_dns cvRoberts_dns.c)
set_target_properties(cvRoberts_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_dns cvRoberts_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_dns.out)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dns.c cvRoberts_dns.out cvRoberts_dns_stats.csv
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_dns)
endif()
