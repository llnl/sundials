# CVODE Kokkos two-dimensional batched Brusselator example
set(examples_list "cv_bruss_batched_kokkos_2D.cpp\;\;develop")
foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})
  set(example_target cv_bruss_batched_kokkos_2D.${backend})
  if(NOT TARGET ${example_target})
    sundials_add_executable(${example_target} cv_bruss_batched_kokkos_2D.cpp)
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")
    target_compile_definitions(${example_target} PRIVATE USE_${backend})
    if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
      target_include_directories(${example_target} PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)
    else()
      target_include_directories(${example_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    target_link_libraries(${example_target}
      PRIVATE sundials_cvode sundials_nveckokkos sundials_sunmatrixkokkosdense
              sundials_sunlinsolkokkosdense ${EXE_EXTRA_LINK_LIBS})
    if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
      target_link_libraries(${example_target} PRIVATE sundials_cvode_fused_stubs)
    endif()
  endif()
  sundials_add_test(${example_target} ${example_target}
                    EXTRA_ARGS --kokkos-disable-warnings
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ${example_target}.out EXAMPLE_TYPE develop)
endforeach()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out cv_bruss_batched_kokkos_2D*.out)
  install(FILES cv_bruss_batched_kokkos_2D.cpp ${example_out}
          "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/cv_bruss_batched_kokkos_2D)
endif()
