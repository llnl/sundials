sundials_add_executable(cvDiurnal_kry_bp cvDiurnal_kry_bp.c)
set_target_properties(cvDiurnal_kry_bp PROPERTIES FOLDER "Examples")
target_link_libraries(cvDiurnal_kry_bp PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvDiurnal_kry_bp cvDiurnal_kry_bp ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvDiurnal_kry_bp.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvDiurnal_kry_bp.c cvDiurnal_kry_bp.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvDiurnal_kry_bp)
endif()
