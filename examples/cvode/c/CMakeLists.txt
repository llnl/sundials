# ---------------------------------------------------------------
# CMakeLists.txt file for CVODE C examples
# ---------------------------------------------------------------

# Each C example has its own source, answer files, and CMakeLists.txt.  Keep
# backend selection here so that the leaf files only describe one executable.
set(CVODE_C_SERIAL_EXAMPLES
    cvAdvDiff_bnd cvAnalytic_mels cvDirectDemo_ls cvDisc_dns cvDiurnal_kry_bp
    cvDiurnal_kry cvKrylovDemo_ls cvKrylovDemo_prec cvParticle_dns
    cvPendulum_dns cvRoberts_dns cvRoberts_dns_constraints cvRoberts_dns_negsol
    cvRoberts_dns_uw cvRocket_dns cvVdp_auto_nls)
set(CVODE_C_LAPACK_EXAMPLES)
set(CVODE_C_KLU_EXAMPLES)
set(CVODE_C_SUPERLUMT_EXAMPLES)
set(CVODE_C_OPENMP_EXAMPLES)
set(CVODE_C_OPENMPDEV_EXAMPLES)
set(CVODE_C_MPI_EXAMPLES)
set(CVODE_C_MPIMANYVECTOR_EXAMPLES)
set(CVODE_C_PARHYP_EXAMPLES)
set(CVODE_C_PETSC_EXAMPLES)
set(CVODE_C_LEAVES ${CVODE_C_SERIAL_EXAMPLES})

set(CVODE_C_LIBS sundials_cvode)
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_C_LIBS sundials_cvode_fused_stubs)
endif()
add_subdirectory(cvAdvDiff_bnd)
add_subdirectory(cvAnalytic_mels)
add_subdirectory(cvDirectDemo_ls)
add_subdirectory(cvDisc_dns)
add_subdirectory(cvDiurnal_kry_bp)
add_subdirectory(cvDiurnal_kry)
add_subdirectory(cvKrylovDemo_ls)
add_subdirectory(cvKrylovDemo_prec)
add_subdirectory(cvParticle_dns)
add_subdirectory(cvPendulum_dns)
add_subdirectory(cvRoberts_dns)
add_subdirectory(cvRoberts_dns_constraints)
add_subdirectory(cvRoberts_dns_negsol)
add_subdirectory(cvRoberts_dns_uw)
add_subdirectory(cvRocket_dns)
add_subdirectory(cvVdp_auto_nls)

if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
   AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
  set(CVODE_C_LAPACK_EXAMPLES cvAdvDiff_bndL cvRoberts_dnsL)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_LAPACK_EXAMPLES})
  add_subdirectory(cvAdvDiff_bndL)
  add_subdirectory(cvRoberts_dnsL)
endif()

if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
  set(CVODE_C_KLU_EXAMPLES cvRoberts_block_klu cvRoberts_klu)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_KLU_EXAMPLES})
  add_subdirectory(cvRoberts_block_klu)
  add_subdirectory(cvRoberts_klu)
endif()

if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
  set(CVODE_C_SUPERLUMT_EXAMPLES cvRoberts_sps)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_SUPERLUMT_EXAMPLES})
  add_subdirectory(cvRoberts_sps)
endif()

if(SUNDIALS_ENABLE_OPENMP)
  set(CVODE_C_OPENMP_EXAMPLES cvAdvDiff_bnd_omp)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_OPENMP_EXAMPLES})
  add_subdirectory(cvAdvDiff_bnd_omp)
endif()

if(SUNDIALS_ENABLE_OPENMP_DEVICE AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
  set(CVODE_C_OPENMPDEV_EXAMPLES cvAdvDiff_kry_ompdev)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_OPENMPDEV_EXAMPLES})
  add_subdirectory(cvAdvDiff_kry_ompdev)
endif()

if(SUNDIALS_ENABLE_MPI)
  set(CVODE_C_MPI_EXAMPLES
      cvAdvDiff_diag_p cvAdvDiff_non_p cvDiurnal_kry_bbd_p cvDiurnal_kry_p)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_MPI_EXAMPLES})
  add_subdirectory(cvAdvDiff_diag_p)
  add_subdirectory(cvAdvDiff_non_p)
  add_subdirectory(cvDiurnal_kry_bbd_p)
  add_subdirectory(cvDiurnal_kry_p)

  if(SUNDIALS_ENABLE_NVECTOR_MPIMANYVECTOR)
    set(CVODE_C_MPIMANYVECTOR_EXAMPLES cvDiurnal_kry_mpimanyvec)
    list(APPEND CVODE_C_LEAVES ${CVODE_C_MPIMANYVECTOR_EXAMPLES})
    add_subdirectory(cvDiurnal_kry_mpimanyvec)
  endif()
endif()

if(SUNDIALS_ENABLE_HYPRE)
  set(CVODE_C_PARHYP_EXAMPLES cvAdvDiff_non_ph)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_PARHYP_EXAMPLES})
  add_subdirectory(cvAdvDiff_non_ph)
endif()

if(SUNDIALS_ENABLE_PETSC)
  set(CVODE_C_PETSC_EXAMPLES cvAdvDiff_petsc cv_petsc_ex7)
  list(APPEND CVODE_C_LEAVES ${CVODE_C_PETSC_EXAMPLES})
  add_subdirectory(cvAdvDiff_petsc)
  add_subdirectory(cv_petsc_ex7)
endif()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES README.md
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c)

  foreach(leaf ${CVODE_C_LEAVES})
    install(FILES ${leaf}/CMakeLists.txt
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c/${leaf})
  endforeach()

  set(CVODE_C_COMPONENTS cvode nvecserial)
  if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
    list(APPEND CVODE_C_COMPONENTS cvode_fused_stubs)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    list(APPEND CVODE_C_COMPONENTS sunlinsollapackband sunlinsollapackdense)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    list(APPEND CVODE_C_COMPONENTS sunlinsolklu)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    list(APPEND CVODE_C_COMPONENTS sunlinsolsuperlumt)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    list(APPEND CVODE_C_COMPONENTS nvecopenmp)
  endif()
  if(SUNDIALS_ENABLE_OPENMP_DEVICE AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    list(APPEND CVODE_C_COMPONENTS nvecopenmpdev)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    list(APPEND CVODE_C_COMPONENTS nvecparallel)
  endif()
  if(SUNDIALS_ENABLE_NVECTOR_MPIMANYVECTOR)
    list(APPEND CVODE_C_COMPONENTS nvecmpimanyvector)
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    list(APPEND CVODE_C_COMPONENTS nvecparhyp)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    list(APPEND CVODE_C_COMPONENTS nvecpetsc sunnonlinsolpetscsnes)
  endif()
  list(REMOVE_DUPLICATES CVODE_C_COMPONENTS)
  list2string(CVODE_C_COMPONENTS CVODE_C_COMPONENTS_STRING)

  if(SUNDIALS_ENABLE_MPI OR SUNDIALS_ENABLE_HYPRE OR SUNDIALS_ENABLE_PETSC)
    set(CVODE_C_NEEDS_MPI ON)
    set(CVODE_C_FIND_MPI "find_package(MPI REQUIRED)")
  else()
    set(CVODE_C_NEEDS_MPI OFF)
    set(CVODE_C_FIND_MPI "")
  endif()

  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_cvode_c_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/c/CMakeLists.txt @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/c/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/c)
  sundials_add_test_install(cvode c EXECUTABLE cvRoberts_dns)
endif()
