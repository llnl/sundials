sundials_add_executable(cvRoberts_dnsL cvRoberts_dnsL.c)
set_target_properties(cvRoberts_dnsL PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_dnsL PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsollapackband sundials_sunlinsollapackdense
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_dnsL cvRoberts_dnsL ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_dnsL.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dnsL.c cvRoberts_dnsL.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_dnsL)
endif()
