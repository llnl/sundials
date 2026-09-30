# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Standalone CMake project shared by installed CVODE example variants.
# -----------------------------------------------------------------

file(
  GLOB
  _cvode_example_sources
  CONFIGURE_DEPENDS
  "*.c"
  "*.cpp"
  "*.cu"
  "*.f90")
foreach(_source IN LISTS _cvode_example_sources)
  if(_source MATCHES "\\.cu$")
    enable_language(CUDA)
  elseif(_source MATCHES "\\.cpp$")
    enable_language(CXX)
  elseif(_source MATCHES "\\.c$")
    enable_language(C)
  elseif(_source MATCHES "\\.f90$")
    enable_language(Fortran)
  endif()
endforeach()

set(CMAKE_Fortran_PREPROCESS ON)
include(CTest)
find_package(SUNDIALS REQUIRED)

if(SUNDIALS_ENABLE_HIP)
  find_package(HIP REQUIRED)
endif()
if((CMAKE_CURRENT_SOURCE_DIR MATCHES "/cpp-magma$" AND SUNDIALS_MAGMA_BACKENDS
                                                       MATCHES "CUDA")
   OR (CMAKE_CURRENT_SOURCE_DIR MATCHES "/cpp-raja$" AND SUNDIALS_RAJA_BACKENDS
                                                         MATCHES "CUDA")
   OR (CMAKE_CURRENT_SOURCE_DIR MATCHES "/cpp-ginkgo$"
       AND SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA"))
  enable_language(CUDA)
elseif(CMAKE_CURRENT_SOURCE_DIR MATCHES "/cpp-magma$"
       AND SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
  enable_language(HIP)
endif()

set(SUNDIALS_ENABLE_EXAMPLES_INSTALL OFF)
set(EXE_EXTRA_LINK_LIBS "${SUNDIALS_MATH_LIBRARY}")
set(SUPERLUMT_LIBRARIES "")

function(sundials_add_executable name)
  set(options)
  set(single_value_args SCALAR_TYPE)
  set(multi_value_args LINK_LIBRARIES INCLUDE_DIRECTORIES COMPILE_DEFINITIONS
                       PROPERTIES)
  cmake_parse_arguments(arg "${options}" "${single_value_args}"
                        "${multi_value_args}" ${ARGN})

  add_executable(${name} ${arg_UNPARSED_ARGUMENTS})
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  if(arg_LINK_LIBRARIES)
    target_link_libraries(${name} ${arg_LINK_LIBRARIES})
  endif()
  if(arg_INCLUDE_DIRECTORIES)
    target_include_directories(${name} ${arg_INCLUDE_DIRECTORIES})
  endif()
  if(arg_COMPILE_DEFINITIONS)
    target_compile_definitions(${name} ${arg_COMPILE_DEFINITIONS})
  endif()
  if(arg_PROPERTIES)
    set_target_properties(${name} PROPERTIES ${arg_PROPERTIES})
  endif()
endfunction()

function(sundials_add_test name executable)
  set(options NODIFF)
  set(one_value_args MPI_NPROCS FLOAT_PRECISION INTEGER_PRECISION ANSWER_DIR
                     ANSWER_FILE EXAMPLE_TYPE)
  set(multi_value_args TEST_ARGS EXTRA_ARGS LABELS)
  cmake_parse_arguments(arg "${options}" "${one_value_args}"
                        "${multi_value_args}" ${ARGN})
  if(arg_MPI_NPROCS)
    add_test(
      NAME ${name}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${arg_MPI_NPROCS}
              $<TARGET_FILE:${executable}> ${arg_TEST_ARGS} ${arg_EXTRA_ARGS})
  else()
    add_test(NAME ${name} COMMAND $<TARGET_FILE:${executable}> ${arg_TEST_ARGS}
                                  ${arg_EXTRA_ARGS})
  endif()
  if(arg_LABELS)
    set_tests_properties(${name} PROPERTIES LABELS "${arg_LABELS}")
  endif()
endfunction()

include(${CMAKE_CURRENT_LIST_DIR}/SundialsAddExample.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/SundialsAddExamplesGinkgo.cmake)

# The build tree and installed package use unnamespaced and namespaced target
# names, respectively. Provide the former locally so the same example files work
# in both contexts.
set(_cvode_targets
    cvode
    cvode_fused_stubs
    cvode_fused_cuda
    cvode_fused_hip
    fcvode_mod
    fnvecserial_mod
    fnvecparallel_mod
    nvecserial
    nvecparallel
    nvecmpimanyvector
    nvecopenmp
    nvecopenmpdev
    nvecparhyp
    nvecpetsc
    nveccuda
    nvechip
    nvecsycl
    nvecraja
    nveckokkos
    sunlinsollapackband
    sunlinsollapackdense
    sunlinsolklu
    sunlinsolsuperlumt
    sunlinsolsuperludist
    sunlinsolcusolversp
    sunnonlinsolpetscsnes
    sunmatrixkokkosdense
    sunlinsolkokkosdense
    sunlinsolonemkldense
    sunlinsolmagmadense
    fsunlinsolklu_mod
    fsunlinsollapackdense_mod)
foreach(_target IN LISTS _cvode_targets)
  if(TARGET SUNDIALS::${_target})
    add_library(sundials_${_target} ALIAS SUNDIALS::${_target})
  endif()
endforeach()

set(CVODE_C_LIBS sundials_cvode)
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_C_LIBS sundials_cvode_fused_stubs)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/cvode-example.cmake)
