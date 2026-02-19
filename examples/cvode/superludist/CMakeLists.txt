# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for CVODE SuperLU_DIST examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;type" where the type is develop
# for examples excluded from 'make test' in releases

list(APPEND CVODE_examples "cvAdvDiff_sludist.cpp\;2 1\;2\;develop")

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  set(_fused_lib sundials_cvode_fused_stubs)
endif()

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # example source files
  sundials_add_executable(${example_target} ${example})

  set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

  if(SUNDIALS_INDEX_SIZE STREQUAL "32")
    set(integer_precision "15")

    # superlu-dist seems to result in more output variability (maybe due to a
    # bug with 32-bit indices?)
    set(float_precision "3")
  else()
    set(integer_precision "default")
    set(float_precision "default")
  endif()

  # add example to regression tests
  sundials_add_test(
    ${example_target} ${example_target}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_nodes}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}_${SUNDIALS_INDEX_SIZE}.out
    EXAMPLE_TYPE ${example_type}
    INTEGER_PRECISION ${integer_precision}
    FLOAT_PRECISION ${float_precision})

  # libraries to link against
  target_link_libraries(
    ${example_target}
    PRIVATE MPI::MPI_CXX sundials_cvode sundials_nvecparallel
            sundials_sunlinsolsuperludist ${_fused_lib} ${EXE_EXTRA_LINK_LIBS})

endforeach(example_tuple ${CVODE_examples})

if(SUPERLUDIST_CUDA)
  set(_ex_lang CUDA)
elseif(SUPERLUDIST_ROCM)
  set(_ex_lang HIP)
else()
  set(_ex_lang CXX)
endif()

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(
    cvode_superludist CVODE_examples
    EXTRA_FILES README
    CMAKE_TEMPLATE cmakelists_${_ex_lang}_MPI_ex.in
    SUNDIALS_TARGETS cvode nvecparallel sunlinsolsuperludist
    DESTINATION cvode/superludist)

  # add test_install target
  sundials_add_test_install(cvode superludist)

endif()
