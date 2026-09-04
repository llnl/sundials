sundials_add_executable(cvRoberts_sps cvRoberts_sps.c)
set_target_properties(cvRoberts_sps PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_sps PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsolsuperlumt ${SUPERLUMT_LIBRARIES}
                      ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  sundials_add_test(cvRoberts_sps cvRoberts_sps ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE cvRoberts_sps.out EXAMPLE_TYPE develop)
endif()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_sps.c cvRoberts_sps.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvRoberts_sps)
endif()
