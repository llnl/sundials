sundials_add_executable(cvAdvDiff_bndL cvAdvDiff_bndL.c)
set_target_properties(cvAdvDiff_bndL PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_bndL PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsollapackband sundials_sunlinsollapackdense
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_bndL cvAdvDiff_bndL ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_bndL.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_bndL.c cvAdvDiff_bndL.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_bndL)
endif()
