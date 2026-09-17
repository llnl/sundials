if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(c)
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran)
endif()
