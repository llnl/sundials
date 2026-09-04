if(MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  include_directories(${MPI_INCLUDE_PATH})
endif()
sundials_add_executable(cvDiurnal_kry_p cvDiurnal_kry_p.c)
set_target_properties(cvDiurnal_kry_p PROPERTIES FOLDER "Examples")
target_link_libraries(cvDiurnal_kry_p PRIVATE ${CVODE_C_LIBS} sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
if(NOT MPI_C_COMPILER)
  target_link_libraries(cvDiurnal_kry_p PRIVATE ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
endif()
sundials_add_test(cvDiurnal_kry_p cvDiurnal_kry_p MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDiurnal_kry_p.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDiurnal_kry_p.c cvDiurnal_kry_p.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvDiurnal_kry_p)
endif()
