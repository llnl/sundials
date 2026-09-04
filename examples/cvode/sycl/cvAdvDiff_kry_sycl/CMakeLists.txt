set(example cvAdvDiff_kry_sycl)
sundials_add_executable(${example} ${example}.cpp)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
target_link_libraries(${example} PRIVATE sundials_cvode sundials_nvecsycl
                          ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(${example} ${example} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ${example}.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES ${example}.cpp ${example}.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/sycl/cvAdvDiff_kry_sycl)
endif()
