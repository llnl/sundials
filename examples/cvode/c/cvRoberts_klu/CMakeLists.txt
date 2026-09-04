sundials_add_executable(cvRoberts_klu cvRoberts_klu.c)
set_target_properties(cvRoberts_klu PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_klu PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_klu cvRoberts_klu ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_klu.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_klu.c cvRoberts_klu.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_klu)
endif()
