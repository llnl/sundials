# ---------------------------------------------------------------
# CMakeLists.txt file for CVODE C++ examples
# ---------------------------------------------------------------

# The C++ examples are organized one executable per directory.  Keep backend
# feature checks here so the leaf files only describe their executable.
add_subdirectory(cv_heat2D)
add_subdirectory(cv_kpr)

if(SUNDIALS_ENABLE_MPI)
  add_subdirectory(cv_heat2D_p)
  if(SUNDIALS_ENABLE_HYPRE)
    add_subdirectory(cv_heat2D_hypre_ls)
    add_subdirectory(cv_heat2D_hypre_pfmg)
  endif()
  if(SUNDIALS_ENABLE_SUPERLUDIST)
    add_subdirectory(cvAdvDiff_sludist)
  endif()
endif()

if(SUNDIALS_ENABLE_MAGMA)
  add_subdirectory(cv_bruss_batched_magma)
endif()

if(SUNDIALS_ENABLE_GINKGO)
  add_subdirectory(cv_bruss_batched_ginkgo)
  add_subdirectory(cv_heat2D_ginkgo)
  add_subdirectory(cv_kpr_ginkgo)
endif()

if(SUNDIALS_ENABLE_SUNLINSOL_KOKKOSDENSE)
  add_subdirectory(cv_bruss_batched_kokkos)
  add_subdirectory(cv_bruss_batched_kokkos_2D)
endif()

if(SUNDIALS_ENABLE_SUNLINSOL_ONEMKLDENSE)
  add_subdirectory(cvRoberts_blockdiag_onemkl)
endif()

if(SUNDIALS_ENABLE_RAJA)
  add_subdirectory(cvAdvDiff_kry_raja)
endif()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  set(CVODE_CPP_COMPONENTS cvode nvecserial)
  set(CVODE_CPP_LEAVES
      cv_heat2D cv_kpr cv_heat2D_p cv_heat2D_hypre_ls cv_heat2D_hypre_pfmg
      cvAdvDiff_sludist cv_bruss_batched_magma cv_bruss_batched_ginkgo
      cv_heat2D_ginkgo cv_kpr_ginkgo cv_bruss_batched_kokkos
      cv_bruss_batched_kokkos_2D cvRoberts_blockdiag_onemkl
      cvAdvDiff_kry_raja)

  if(SUNDIALS_ENABLE_MPI)
    list(APPEND CVODE_CPP_COMPONENTS nvecparallel)
    set(CVODE_CPP_FIND_MPI "find_package(MPI REQUIRED)")
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    list(APPEND CVODE_CPP_COMPONENTS nvecparallel)
  endif()
  if(SUNDIALS_ENABLE_SUPERLUDIST)
    list(APPEND CVODE_CPP_COMPONENTS nvecparallel sunlinsolsuperludist)
  endif()
  if(SUNDIALS_ENABLE_MAGMA)
    list(APPEND CVODE_CPP_COMPONENTS sunlinsolmagmadense)
    if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA")
      list(APPEND CVODE_CPP_COMPONENTS nveccuda)
    elseif(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
      list(APPEND CVODE_CPP_COMPONENTS nvechip)
    endif()
  endif()
  if(SUNDIALS_ENABLE_GINKGO)
    if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
      list(APPEND CVODE_CPP_COMPONENTS nveccuda)
    endif()
    if(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
      list(APPEND CVODE_CPP_COMPONENTS nvechip)
    endif()
    if(SUNDIALS_GINKGO_BACKENDS MATCHES "SYCL")
      list(APPEND CVODE_CPP_COMPONENTS nvecsycl)
    endif()
    if(SUNDIALS_GINKGO_BACKENDS MATCHES "OMP")
      list(APPEND CVODE_CPP_COMPONENTS nvecopenmp)
    endif()
    if(SUNDIALS_GINKGO_BACKENDS MATCHES "REF")
      list(APPEND CVODE_CPP_COMPONENTS nvecserial)
    endif()
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KOKKOSDENSE)
    list(APPEND CVODE_CPP_COMPONENTS nveckokkos sunmatrixkokkosdense
         sunlinsolkokkosdense)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_ONEMKLDENSE)
    list(APPEND CVODE_CPP_COMPONENTS nvecsycl sunlinsolonemkldense)
  endif()
  if(SUNDIALS_ENABLE_RAJA)
    list(APPEND CVODE_CPP_COMPONENTS nvecraja)
  endif()
  if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
    list(APPEND CVODE_CPP_COMPONENTS cvode_fused_stubs)
  endif()
  list(REMOVE_DUPLICATES CVODE_CPP_COMPONENTS)

  set(CVODE_CPP_LANGUAGES)
  if((SUNDIALS_ENABLE_MAGMA AND SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA") OR
     (SUNDIALS_ENABLE_RAJA AND SUNDIALS_RAJA_BACKENDS MATCHES "CUDA") OR
     (SUNDIALS_ENABLE_GINKGO AND SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA") OR
     (SUNDIALS_ENABLE_SUPERLUDIST AND SUPERLUDIST_CUDA))
    set(CVODE_CPP_LANGUAGES CUDA)
  endif()

  install(FILES README.md
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp)

  foreach(leaf ${CVODE_CPP_LEAVES})
    install(FILES ${leaf}/CMakeLists.txt
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp/${leaf})
  endforeach()

  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_cvode_cpp_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/cpp/CMakeLists.txt @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/cpp/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cpp)
  sundials_add_test_install(cvode cpp)
endif()
