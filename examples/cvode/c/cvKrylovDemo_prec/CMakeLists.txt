sundials_add_executable(cvKrylovDemo_prec cvKrylovDemo_prec.c)
set_target_properties(cvKrylovDemo_prec PROPERTIES FOLDER "Examples")
target_link_libraries(cvKrylovDemo_prec PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvKrylovDemo_prec cvKrylovDemo_prec ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvKrylovDemo_prec.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvKrylovDemo_prec.c cvKrylovDemo_prec.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvKrylovDemo_prec)
endif()
