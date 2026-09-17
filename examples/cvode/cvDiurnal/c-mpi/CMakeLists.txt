sundials_add_executable(cvDiurnal_kry_bbd_p cvDiurnal_kry_bbd_p.c)
set_target_properties(cvDiurnal_kry_bbd_p PROPERTIES FOLDER "Examples")
target_link_libraries(cvDiurnal_kry_bbd_p PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDiurnal_kry_bbd_p cvDiurnal_kry_bbd_p MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDiurnal_kry_bbd_p.out EXAMPLE_TYPE exclude-sigle)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDiurnal_kry_bbd_p.c cvDiurnal_kry_bbd_p.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvDiurnal/c-mpi)
endif()

sundials_add_executable(cvDiurnal_kry_p cvDiurnal_kry_p.c)
set_target_properties(cvDiurnal_kry_p PROPERTIES FOLDER "Examples")
target_link_libraries(cvDiurnal_kry_p PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDiurnal_kry_p cvDiurnal_kry_p MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDiurnal_kry_p.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDiurnal_kry_p.c cvDiurnal_kry_p.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvDiurnal/c-mpi)
endif()
