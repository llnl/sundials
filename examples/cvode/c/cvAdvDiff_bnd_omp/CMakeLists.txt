sundials_add_executable(cvAdvDiff_bnd_omp cvAdvDiff_bnd_omp.c)
set_target_properties(cvAdvDiff_bnd_omp PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_bnd_omp PRIVATE ${CVODE_C_LIBS} sundials_nvecopenmp
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_bnd_omp cvAdvDiff_bnd_omp TEST_ARGS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_bnd_omp.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_bnd_omp.c cvAdvDiff_bnd_omp.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_bnd_omp)
endif()
