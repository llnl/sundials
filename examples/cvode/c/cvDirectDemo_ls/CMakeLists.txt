sundials_add_executable(cvDirectDemo_ls cvDirectDemo_ls.c)
set_target_properties(cvDirectDemo_ls PROPERTIES FOLDER "Examples")
target_link_libraries(cvDirectDemo_ls PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDirectDemo_ls cvDirectDemo_ls ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDirectDemo_ls.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDirectDemo_ls.c cvDirectDemo_ls.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvDirectDemo_ls)
endif()
