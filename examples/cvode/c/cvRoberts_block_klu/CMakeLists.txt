sundials_add_executable(cvRoberts_block_klu cvRoberts_block_klu.c)
set_target_properties(cvRoberts_block_klu PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_block_klu PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_block_klu cvRoberts_block_klu
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_block_klu.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_block_klu.c cvRoberts_block_klu.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_block_klu)
endif()
