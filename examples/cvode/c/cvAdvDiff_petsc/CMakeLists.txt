if(MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  include_directories(${MPI_INCLUDE_PATH})
endif()
sundials_add_executable(cvAdvDiff_petsc cvAdvDiff_petsc.c)
set_target_properties(cvAdvDiff_petsc PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_petsc PRIVATE ${CVODE_C_LIBS} sundials_nvecpetsc
                      sundials_sunnonlinsolpetscsnes ${EXE_EXTRA_LINK_LIBS})
if(NOT MPI_C_COMPILER)
  target_link_libraries(cvAdvDiff_petsc PRIVATE ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
endif()
sundials_add_test(cvAdvDiff_petsc cvAdvDiff_petsc MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_petsc.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_petsc.c cvAdvDiff_petsc.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_petsc)
endif()
