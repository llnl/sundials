sundials_add_executable(cvAdvDiff_bnd cvAdvDiff_bnd.c)
set_target_properties(cvAdvDiff_bnd PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_bnd PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_bnd cvAdvDiff_bnd ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_bnd.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_bnd.c cvAdvDiff_bnd.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_bnd)
endif()
