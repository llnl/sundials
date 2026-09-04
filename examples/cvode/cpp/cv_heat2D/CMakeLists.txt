# CVODE serial C++ heat equation example
set(CVODE_examples "cv_heat2D.cpp\;\;\;\;exclude-single")

sundials_add_executable(cv_heat2D cv_heat2D.cpp)
set_target_properties(cv_heat2D PROPERTIES FOLDER "Examples")
if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
  target_include_directories(cv_heat2D PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)
else()
  target_include_directories(cv_heat2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
target_link_libraries(cv_heat2D PRIVATE sundials_cvode sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_heat2D PRIVATE sundials_cvode_fused_stubs)
endif()
sundials_add_test(cv_heat2D cv_heat2D ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_heat2D.out EXAMPLE_TYPE exclude-single)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_heat2D.cpp cv_heat2D.hpp cv_heat2D.out plot_heat2D.py
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cv_heat2D)
endif()
