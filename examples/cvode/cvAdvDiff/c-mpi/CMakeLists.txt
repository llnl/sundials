sundials_add_executable(cvAdvDiff_diag_p cvAdvDiff_diag_p.c)
set_target_properties(cvAdvDiff_diag_p PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_diag_p PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_diag_p cvAdvDiff_diag_p MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_diag_p.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_diag_p.c cvAdvDiff_diag_p.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvAdvDiff/c-mpi)
endif()

sundials_add_executable(cvAdvDiff_non_p cvAdvDiff_non_p.c)
set_target_properties(cvAdvDiff_non_p PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_non_p PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_non_p cvAdvDiff_non_p MPI_NPROCS 2
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_non_p.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_non_p.c cvAdvDiff_non_p.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvAdvDiff/c-mpi)
endif()
