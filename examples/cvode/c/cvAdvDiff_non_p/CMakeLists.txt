if(MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  include_directories(${MPI_INCLUDE_PATH})
endif()
sundials_add_executable(cvAdvDiff_non_p cvAdvDiff_non_p.c)
set_target_properties(cvAdvDiff_non_p PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_non_p PRIVATE ${CVODE_C_LIBS} sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
if(NOT MPI_C_COMPILER)
  target_link_libraries(cvAdvDiff_non_p PRIVATE ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
endif()
sundials_add_test(cvAdvDiff_non_p cvAdvDiff_non_p MPI_NPROCS 2
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_non_p.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_non_p.c cvAdvDiff_non_p.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvAdvDiff_non_p)
endif()
