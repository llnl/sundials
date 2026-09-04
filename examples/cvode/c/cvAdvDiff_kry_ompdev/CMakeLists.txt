sundials_add_executable(cvAdvDiff_kry_ompdev cvAdvDiff_kry_ompdev.c)
set_target_properties(cvAdvDiff_kry_ompdev PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_kry_ompdev PRIVATE ${CVODE_C_LIBS}
                      sundials_nvecopenmpdev ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_kry_ompdev cvAdvDiff_kry_ompdev TEST_ARGS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_kry_ompdev.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_kry_ompdev.c cvAdvDiff_kry_ompdev.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_kry_ompdev)
endif()
