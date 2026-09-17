sundials_add_executable(cvDiurnal_kry_mpimanyvec cvDiurnal_kry_mpimanyvec.c)
set_target_properties(cvDiurnal_kry_mpimanyvec PROPERTIES FOLDER "Examples")
target_link_libraries(cvDiurnal_kry_mpimanyvec PRIVATE MPI::MPI_C ${CVODE_C_LIBS}
                      sundials_nvecmpimanyvector sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDiurnal_kry_mpimanyvec cvDiurnal_kry_mpimanyvec MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDiurnal_kry_mpimanyvec.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDiurnal_kry_mpimanyvec.c cvDiurnal_kry_mpimanyvec.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvDiurnal/c-mpi-mpimanyvector)
endif()
