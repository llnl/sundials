# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODE Fortran examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases.
set(FCVODE_examples
    "cv_advdiff_bnd_f2003\;develop")
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  list(APPEND FCVODE_examples
       "cv_analytic_fp_f2003\;develop"
       "cv_analytic_sys_dns_f2003\;develop"
       "cv_analytic_sys_dns_jac_f2003\;develop"
       "cv_brusselator_dns_f2003\;develop"
       "cv_diurnal_kry_f2003\;develop"
       "cv_diurnal_kry_bp_f2003\;develop"
       "cv_roberts_dns_constraints_f2003\;develop"
       "cv_roberts_dns_f2003\;develop")
  set(FCVODE_examples_KLU
      "cv_analytic_sys_klu_f2003\;develop"
      "cv_roberts_klu_f2003\;develop")
  set(FCVODE_examples_LAPACK "cv_roberts_dnsL_f2003\;develop")
else()
  set(FCVODE_examples_KLU "")
  set(FCVODE_examples_LAPACK "")
endif()

# Parallel examples are available with 64-bit indices and MPI.
set(FCVODE_examples_parallel
    "cv_diag_kry_bbd_f2003\;\;1\;4\;develop\;2"
    "cv_diag_kry_f2003\;\;1\;4\;develop\;2"
    "cv_diag_non_p_f2003\;\;1\;4\;develop\;2")

add_subdirectory(cv_advdiff_bnd_f2003)
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(cv_analytic_fp_f2003)
  add_subdirectory(cv_analytic_sys_dns_f2003)
  add_subdirectory(cv_analytic_sys_dns_jac_f2003)
  add_subdirectory(cv_brusselator_dns_f2003)
  add_subdirectory(cv_diurnal_kry_f2003)
  add_subdirectory(cv_diurnal_kry_bp_f2003)
  add_subdirectory(cv_roberts_dns_constraints_f2003)
  add_subdirectory(cv_roberts_dns_f2003)
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(cv_analytic_sys_klu_f2003)
    add_subdirectory(cv_roberts_klu_f2003)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(cv_roberts_dnsL_f2003)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(cv_diag_kry_bbd_f2003)
    add_subdirectory(cv_diag_kry_f2003)
    add_subdirectory(cv_diag_non_p_f2003)
  endif()
endif()

# Keep the installed source hierarchy aligned with the source tree, while
# placing all Fortran examples under cvode/fortran.
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  set(SOLVER CVODE)
  set(SOLVER_LIB sundials_cvode)
  set(SOLVER_FLIB sundials_fcvode_mod)
  set(NVEC_LIB sundials_nvecserial)
  set(NVEC_FLIB sundials_fnvecserial_mod)
  set(FC_COMPILER ${_EXAMPLES_Fortran_COMPILER})
  set(MPI_LIBRARIES "")
  if(SUNDIALS_ENABLE_MPI)
    set(FC_COMPILER ${MPI_Fortran_COMPILER})
    list2string(MPI_Fortran_LIBRARIES MPI_LIBRARIES)
  endif()
  set(LIBS "")
  if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs")
  endif()
  set(FCVODE_examples_install
      "cv_advdiff_bnd_f2003/cv_advdiff_bnd_f2003\;develop")
  set(FCVODE_examples_klu_install "")
  set(FCVODE_examples_lapack_install "")
  if(SUNDIALS_INDEX_SIZE MATCHES "64")
    list(APPEND FCVODE_examples_install
         "cv_analytic_fp_f2003/cv_analytic_fp_f2003\;develop"
         "cv_analytic_sys_dns_f2003/cv_analytic_sys_dns_f2003\;develop"
         "cv_analytic_sys_dns_jac_f2003/cv_analytic_sys_dns_jac_f2003\;develop"
         "cv_brusselator_dns_f2003/cv_brusselator_dns_f2003\;develop"
         "cv_diurnal_kry_f2003/cv_diurnal_kry_f2003\;develop"
         "cv_diurnal_kry_bp_f2003/cv_diurnal_kry_bp_f2003\;develop"
         "cv_roberts_dns_constraints_f2003/cv_roberts_dns_constraints_f2003\;develop"
         "cv_roberts_dns_f2003/cv_roberts_dns_f2003\;develop")
    if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
      list(APPEND FCVODE_examples_klu_install
           "cv_analytic_sys_klu_f2003/cv_analytic_sys_klu_f2003\;develop"
           "cv_roberts_klu_f2003/cv_roberts_klu_f2003\;develop")
    endif()
    if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
       AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
      list(APPEND FCVODE_examples_lapack_install
           "cv_roberts_dnsL_f2003/cv_roberts_dnsL_f2003\;develop")
    endif()
  endif()
  examples2string(FCVODE_examples_install EXAMPLES)
  examples2string(FCVODE_examples_klu_install EXAMPLES_KLU)
  examples2string(FCVODE_examples_lapack_install EXAMPLES_LAPACK)
  set(PARALLEL_EXAMPLES "")
  if(SUNDIALS_INDEX_SIZE MATCHES "64" AND SUNDIALS_ENABLE_MPI)
    set(FCVODE_parallel_install
        "cv_diag_kry_bbd_f2003/cv_diag_kry_bbd_f2003\;develop"
        "cv_diag_kry_f2003/cv_diag_kry_f2003\;develop"
        "cv_diag_non_p_f2003/cv_diag_non_p_f2003\;develop")
    examples2string(FCVODE_parallel_install PARALLEL_EXAMPLES)
  endif()
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_nested_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/fortran/CMakeLists.txt @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/fortran/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/fortran)
  if(UNIX)
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_nested_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fortran/Makefile_ex @ONLY)
    install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/fortran/Makefile_ex
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/fortran
            RENAME Makefile)
  endif()
  sundials_add_test_install(cvode fortran)
endif()
