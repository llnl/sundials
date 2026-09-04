set(example cv_diag_kry_f2003)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)
sundials_add_executable(${example} ${example}.f90)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
target_link_libraries(${example} MPI::MPI_Fortran sundials_cvode sundials_fcvode_mod
                       sundials_nvecparallel sundials_fnvecparallel_mod
                       ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(${example} ${example} MPI_NPROCS 4 ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ${example}.out EXAMPLE_TYPE develop FLOAT_PRECISION 2)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES ${example}.f90 ${example}.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/fortran/cv_diag_kry_f2003)
endif()
