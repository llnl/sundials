if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(c)
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(c-lapack)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c-klu)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(c-superlumt)
  endif()
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_SUNLINSOL_ONEMKLDENSE)
  add_subdirectory(cpp-onemkl)
endif()

if(SUNDIALS_ENABLE_CUDA_EXAMPLES AND SUNDIALS_ENABLE_CUDA
   AND CMAKE_CUDA_COMPILER AND SUNDIALS_INDEX_SIZE MATCHES "32")
  add_subdirectory(cuda-cusolversp)
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran)
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(fortran-lapack)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(fortran-klu)
  endif()
endif()
