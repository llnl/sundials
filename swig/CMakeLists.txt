# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# CMake build for generating the SUNDIALS SWIG-Fortran interfaces
# ---------------------------------------------------------------

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.18)
  project(SUNDIALS_SWIG NONE)
  set(SUNDIALS_SWIG_STANDALONE TRUE)
else()
  set(SUNDIALS_SWIG_STANDALONE FALSE)
endif()

if(SUNDIALS_SWIG_STANDALONE
   AND ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}"))
  message(FATAL_ERROR "In-source build prohibited.")
endif()

find_program(SWIG_EXECUTABLE NAMES swig REQUIRED)
message(STATUS "Using SWIG executable: ${SWIG_EXECUTABLE}")

get_filename_component(SUNDIALS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

set(SUNDIALS_SWIG_INCLUDE_DIR
    "${SUNDIALS_ROOT}/include"
    CACHE PATH "SUNDIALS include directory passed to SWIG")

# this list is "<path to swig file>|output path|output name"
set(SUNDIALS_SWIG_INTERFACE_SPECS
    "sundials/fsundials_core_mod.i|src/sundials|fsundials_core_mod"
    "arkode/farkode_mod.i|src/arkode|farkode_mod"
    "arkode/farkode_arkstep_mod.i|src/arkode|farkode_arkstep_mod"
    "arkode/farkode_erkstep_mod.i|src/arkode|farkode_erkstep_mod"
    "arkode/farkode_sprkstep_mod.i|src/arkode|farkode_sprkstep_mod"
    "arkode/farkode_mristep_mod.i|src/arkode|farkode_mristep_mod"
    "arkode/farkode_lsrkstep_mod.i|src/arkode|farkode_lsrkstep_mod"
    "arkode/farkode_splittingstep_mod.i|src/arkode|farkode_splittingstep_mod"
    "arkode/farkode_forcingstep_mod.i|src/arkode|farkode_forcingstep_mod"
    "cvode/fcvode_mod.i|src/cvode|fcvode_mod"
    "cvodes/fcvodes_mod.i|src/cvodes|fcvodes_mod"
    "ida/fida_mod.i|src/ida|fida_mod"
    "idas/fidas_mod.i|src/idas|fidas_mod"
    "kinsol/fkinsol_mod.i|src/kinsol|fkinsol_mod"
    "nvector/fnvector_serial_mod.i|src/nvector/serial|fnvector_serial_mod"
    "nvector/fnvector_parallel_mod.i|src/nvector/parallel|fnvector_parallel_mod"
    "nvector/fnvector_openmp_mod.i|src/nvector/openmp|fnvector_openmp_mod"
    "nvector/fnvector_pthreads_mod.i|src/nvector/pthreads|fnvector_pthreads_mod"
    "nvector/fnvector_manyvector_mod.i|src/nvector/manyvector|fnvector_manyvector_mod"
    "nvector/fnvector_mpimanyvector_mod.i|src/nvector/manyvector|fnvector_mpimanyvector_mod"
    "nvector/fnvector_mpiplusx_mod.i|src/nvector/mpiplusx|fnvector_mpiplusx_mod"
    "sunmatrix/fsunmatrix_band_mod.i|src/sunmatrix/band|fsunmatrix_band_mod"
    "sunmatrix/fsunmatrix_dense_mod.i|src/sunmatrix/dense|fsunmatrix_dense_mod"
    "sunmatrix/fsunmatrix_sparse_mod.i|src/sunmatrix/sparse|fsunmatrix_sparse_mod"
    "sunlinsol/fsunlinsol_band_mod.i|src/sunlinsol/band|fsunlinsol_band_mod"
    "sunlinsol/fsunlinsol_dense_mod.i|src/sunlinsol/dense|fsunlinsol_dense_mod"
    "sunlinsol/fsunlinsol_lapackdense_mod.i|src/sunlinsol/lapackdense|fsunlinsol_lapackdense_mod"
    "sunlinsol/fsunlinsol_klu_mod.i|src/sunlinsol/klu|fsunlinsol_klu_mod"
    "sunlinsol/fsunlinsol_spbcgs_mod.i|src/sunlinsol/spbcgs|fsunlinsol_spbcgs_mod"
    "sunlinsol/fsunlinsol_spfgmr_mod.i|src/sunlinsol/spfgmr|fsunlinsol_spfgmr_mod"
    "sunlinsol/fsunlinsol_spgmr_mod.i|src/sunlinsol/spgmr|fsunlinsol_spgmr_mod"
    "sunlinsol/fsunlinsol_sptfqmr_mod.i|src/sunlinsol/sptfqmr|fsunlinsol_sptfqmr_mod"
    "sunlinsol/fsunlinsol_pcg_mod.i|src/sunlinsol/pcg|fsunlinsol_pcg_mod"
    "sunnonlinsol/fsunnonlinsol_newton_mod.i|src/sunnonlinsol/newton|fsunnonlinsol_newton_mod"
    "sunnonlinsol/fsunnonlinsol_fixedpoint_mod.i|src/sunnonlinsol/fixedpoint|fsunnonlinsol_fixedpoint_mod"
    "sunadaptcontroller/fsunadaptcontroller_imexgus_mod.i|src/sunadaptcontroller/imexgus|fsunadaptcontroller_imexgus_mod"
    "sunadaptcontroller/fsunadaptcontroller_soderlind_mod.i|src/sunadaptcontroller/soderlind|fsunadaptcontroller_soderlind_mod"
    "sunadaptcontroller/fsunadaptcontroller_mrihtol_mod.i|src/sunadaptcontroller/mrihtol|fsunadaptcontroller_mrihtol_mod"
    "sunadjointcheckpointscheme/fsunadjointcheckpointscheme_fixed_mod.i|src/sunadjointcheckpointscheme/fixed|fsunadjointcheckpointscheme_fixed_mod"
    "sundomeigest/fsundomeigest_power_mod.i|src/sundomeigest/power|fsundomeigest_power_mod"
    "sundomeigest/fsundomeigest_arnoldi_mod.i|src/sundomeigest/arnoldi|fsundomeigest_arnoldi_mod")

function(sundials_configure_swig_types int_size scalar_type precision out_var)
  set(PRECISION_LEVEL "#define SUNDIALS_${precision}_PRECISION")
  set(SCALAR_TYPE "#define SUNDIALS_SCALAR_TYPE_${scalar_type}")
  set(INDEX_TYPE "#define SUNDIALS_INT${int_size}_T")

  if("${int_size}" STREQUAL "32")
    set(SUNDIALS_CINDEX_TYPE "int32_t")
  else()
    set(SUNDIALS_CINDEX_TYPE "int64_t")
  endif()

  set(SUNDIALS_COUNTER_TYPE "long int")

  set(_configured_swig_dir "${CMAKE_CURRENT_BINARY_DIR}/sundials")
  file(MAKE_DIRECTORY "${_configured_swig_dir}")

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sundials/fsundials_types.in"
                 "${_configured_swig_dir}/fsundials_types.i" @ONLY)

  set(${out_var}
      "${_configured_swig_dir}/fsundials_types.i"
      PARENT_SCOPE)
endfunction()

function(sundials_collect_swig_outputs int_size scalar_type swig_types_file out_var)
  set(generated_files)

  foreach(spec IN LISTS SUNDIALS_SWIG_INTERFACE_SPECS)
    string(REPLACE "|" ";" fields "${spec}")
    list(GET fields 0 interface_rel)
    list(GET fields 1 output_rel)
    list(GET fields 2 output_name)

    set(interface_file "${CMAKE_CURRENT_SOURCE_DIR}/${interface_rel}")
    set(output_dir "${SUNDIALS_ROOT}/${output_rel}/fmod_int${int_size}")
    set(output_c "${output_dir}/${output_name}.c")
    set(output_f90 "${output_dir}/${output_name}.f90")

    add_custom_command(
      OUTPUT ${output_c} ${output_f90}
      COMMAND
        ${SWIG_EXECUTABLE} -fortran -outdir ${output_dir} -o
        ${output_c} -I${SUNDIALS_SWIG_INCLUDE_DIR}
        -I${CMAKE_CURRENT_BINARY_DIR}/sundials ${interface_file}
      MAIN_DEPENDENCY ${interface_file}
      DEPENDS ${swig_types_file}
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      COMMENT
        "Generating ${output_name} with SWIG"
      VERBATIM)

    list(APPEND generated_files ${output_c} ${output_f90})
  endforeach()

  set(${out_var}
      "${generated_files}"
      PARENT_SCOPE)
endfunction()

if(SUNDIALS_SWIG_STANDALONE)
  set(
    SUNDIALS_INDEX_SIZE
    "64"
    CACHE STRING "SUNDIALS index size to use when regenerating SWIG interfaces")
  set_property(CACHE SUNDIALS_INDEX_SIZE PROPERTY STRINGS 32 64)

  set(
    SUNDIALS_SCALAR_TYPE
    "REAL"
    CACHE STRING
          "SUNDIALS scalar type to use when regenerating SWIG interfaces")
  set_property(CACHE SUNDIALS_SCALAR_TYPE PROPERTY STRINGS REAL COMPLEX)

  set(
    SUNDIALS_PRECISION
    "DOUBLE"
    CACHE STRING "SUNDIALS precision to use when regenerating SWIG interfaces")
  set_property(CACHE SUNDIALS_PRECISION PROPERTY STRINGS SINGLE DOUBLE EXTENDED)

  message(STATUS "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}")
  message(STATUS "SUNDIALS_SCALAR_TYPE=${SUNDIALS_SCALAR_TYPE}")
  message(STATUS "SUNDIALS_PRECISION=${SUNDIALS_PRECISION}")

  sundials_configure_swig_types("${SUNDIALS_INDEX_SIZE}"
                                "${SUNDIALS_SCALAR_TYPE}"
                                "${SUNDIALS_PRECISION}"
                                _sundials_swig_types_file)

  sundials_collect_swig_outputs("${SUNDIALS_INDEX_SIZE}"
                                "${SUNDIALS_SCALAR_TYPE}"
                                "${_sundials_swig_types_file}"
                                _sundials_swig_outputs)
  add_custom_target(swig ALL DEPENDS ${_sundials_swig_outputs})

  message(
    STATUS
      "Added target 'swig'"
  )
else()
  if(NOT DEFINED SUNDIALS_INDEX_SIZE)
    message(FATAL_ERROR "SUNDIALS_INDEX_SIZE must be set when including swig/.")
  endif()

  if(NOT DEFINED SUNDIALS_SCALAR_TYPE)
    message(FATAL_ERROR
            "SUNDIALS_SCALAR_TYPE must be set when including swig/.")
  endif()

  if(NOT DEFINED SUNDIALS_PRECISION)
    message(FATAL_ERROR "SUNDIALS_PRECISION must be set when including swig/.")
  endif()

  sundials_configure_swig_types("${SUNDIALS_INDEX_SIZE}"
                                "${SUNDIALS_SCALAR_TYPE}"
                                "${SUNDIALS_PRECISION}"
                                _sundials_swig_types_file)

  sundials_collect_swig_outputs("${SUNDIALS_INDEX_SIZE}"
                                "${SUNDIALS_SCALAR_TYPE}"
                                "${_sundials_swig_types_file}"
                                _sundials_swig_outputs)
  add_custom_target(sundials_swig ALL DEPENDS ${_sundials_swig_outputs})

  message(
    STATUS
      "Added target 'sundials_swig'"
  )
endif()
