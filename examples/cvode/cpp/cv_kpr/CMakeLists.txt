# CVODE serial C++ KPR example
set(CVODE_examples
    "cv_kpr.cpp\;\;\;\;develop"
    "cv_kpr.cpp\;cvode.init_step 0.01\;\;\;develop")

sundials_add_executable(cv_kpr cv_kpr.cpp)
set_target_properties(cv_kpr PROPERTIES FOLDER "Examples")
if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
  target_include_directories(cv_kpr PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)
else()
  target_include_directories(cv_kpr PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
target_link_libraries(cv_kpr PRIVATE sundials_cvode sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_kpr PRIVATE sundials_cvode_fused_stubs)
endif()
foreach(_tuple ${CVODE_examples})
  list(GET _tuple 1 _args)
  list(GET _tuple 4 _type)
  if(_args STREQUAL "")
    set(_test cv_kpr)
  else()
    string(REPLACE " " "_" _test "cv_kpr_${_args}")
  endif()
  sundials_add_test(${_test} cv_kpr TEST_ARGS ${_args}
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ${_test}.out EXAMPLE_TYPE ${_type})
endforeach()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_kpr.cpp cv_kpr.hpp cv_kpr.out cv_kpr_cvode.init_step_0.01.out
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cv_kpr)
endif()
