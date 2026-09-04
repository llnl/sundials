sundials_add_executable(cvAdvDiff_non_ph cvAdvDiff_non_ph.c)
set_target_properties(cvAdvDiff_non_ph PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_non_ph PRIVATE ${CVODE_C_LIBS} sundials_nvecparhyp
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_non_ph cvAdvDiff_non_ph MPI_NPROCS 2
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_non_ph.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_non_ph.c cvAdvDiff_non_ph.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_non_ph)
endif()
