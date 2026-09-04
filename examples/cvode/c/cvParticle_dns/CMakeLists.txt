sundials_add_executable(cvParticle_dns cvParticle_dns.c)
set_target_properties(cvParticle_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cvParticle_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvParticle_dns cvParticle_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvParticle_dns.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvParticle_dns.c cvParticle_dns.out plot_cvParticle.py
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/cvParticle_dns)
endif()
