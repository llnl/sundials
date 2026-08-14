# -----------------------------------------------------------------------------
# Programmer(s): Steven Smith and Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Module to find and setup KLU.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# KLU does not support single or extended precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR "KLU is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(KLU REQUIRED)

message(STATUS "KLU_LIBRARIES:   ${KLU_LIBRARIES}")
message(STATUS "KLU_INCLUDE_DIR: ${KLU_INCLUDE_DIR}")

set(_KLU_INCLUDE_DIRS)
if(KLU_INCLUDE_DIR)
  list(APPEND _KLU_INCLUDE_DIRS "${KLU_INCLUDE_DIR}")
endif()

if(TARGET SUNDIALS::KLU)
  set(_KLU_TARGET SUNDIALS::KLU)
  get_target_property(_KLU_ALIASED_TARGET SUNDIALS::KLU ALIASED_TARGET)
  if(_KLU_ALIASED_TARGET)
    set(_KLU_TARGET "${_KLU_ALIASED_TARGET}")
  endif()

  foreach(_KLU_INCLUDE_PROPERTY
          INTERFACE_INCLUDE_DIRECTORIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    get_target_property(_KLU_TARGET_INCLUDE_DIRS "${_KLU_TARGET}"
                        "${_KLU_INCLUDE_PROPERTY}")
    if(_KLU_TARGET_INCLUDE_DIRS)
      foreach(_KLU_TARGET_INCLUDE_DIR ${_KLU_TARGET_INCLUDE_DIRS})
        if(NOT _KLU_TARGET_INCLUDE_DIR MATCHES "^\\$<")
          list(APPEND _KLU_INCLUDE_DIRS "${_KLU_TARGET_INCLUDE_DIR}")
        endif()
      endforeach()
    endif()
  endforeach()
endif()

if(_KLU_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES _KLU_INCLUDE_DIRS)
endif()
message(STATUS "KLU_INCLUDE_DIRS: ${_KLU_INCLUDE_DIRS}")

set(_KLU_LINK_LIBRARIES)
if(TARGET SUNDIALS::KLU)
  if(_KLU_TARGET)
    list(APPEND _KLU_LINK_LIBRARIES "${_KLU_TARGET}")
  else()
    list(APPEND _KLU_LINK_LIBRARIES SUNDIALS::KLU)
  endif()
elseif(KLU_LIBRARIES)
  list(APPEND _KLU_LINK_LIBRARIES ${KLU_LIBRARIES})
endif()

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_KLU_CHECKS)

  message(CHECK_START "Testing KLU")

  if(NOT _KLU_INCLUDE_DIRS)
    message(CHECK_FAIL "failed")
    message(
      FATAL_ERROR
        "Could not determine KLU include directories from KLU_INCLUDE_DIR or SUNDIALS::KLU"
    )
  endif()
  if(NOT _KLU_LINK_LIBRARIES)
    message(CHECK_FAIL "failed")
    message(
      FATAL_ERROR
        "Could not determine KLU link libraries from KLU_LIBRARIES or SUNDIALS::KLU"
    )
  endif()

  if(SUNDIALS_INDEX_SIZE MATCHES "64")
    # Check size of SuiteSparse_long
    include(CheckTypeSize)
    set(save_CMAKE_EXTRA_INCLUDE_FILES ${CMAKE_EXTRA_INCLUDE_FILES})
    list(APPEND CMAKE_EXTRA_INCLUDE_FILES "klu.h")
    set(save_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    list(APPEND CMAKE_REQUIRED_INCLUDES ${_KLU_INCLUDE_DIRS})
    check_type_size("SuiteSparse_long" SIZEOF_SUITESPARSE_LONG)
    set(CMAKE_EXTRA_INCLUDE_FILES ${save_CMAKE_EXTRA_INCLUDE_FILES})
    set(CMAKE_REQUIRED_INCLUDES ${save_CMAKE_REQUIRED_INCLUDES})
    message(STATUS "Size of SuiteSparse_long is ${SIZEOF_SUITESPARSE_LONG}")
    if(NOT SIZEOF_SUITESPARSE_LONG EQUAL "8")
      message(CHECK_FAIL "failed")
      message(
        FATAL_ERROR
          "Size of 'sunindextype' is 8 but size of 'SuiteSparse_long' is ${SIZEOF_SUITESPARSE_LONG}. KLU cannot be used."
      )
    endif()
  endif()

  # Create the test directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/KLU_TEST)

  # Create a C source file which calls a KLU function
  file(WRITE ${TEST_DIR}/test.c
       "#include <klu.h>\n"
       "int main(void) {\n"
       "  klu_common Common; (void)Common;\n"
       "  return klu_defaults(&Common) ? 0 : 1;\n"
       "}\n")

  # Link against the resolved imported target when available. Avoid passing a
  # SUNDIALS::KLU alias to try_compile because alias targets to imported
  # libraries require newer CMake versions than SUNDIALS currently requires.
  string(REPLACE ";" "\\;" _KLU_INCLUDE_DIRS_ARG "${_KLU_INCLUDE_DIRS}")
  set(TEST_BINARY_DIR "${TEST_DIR}/build")
  try_compile(
    COMPILE_OK
    "${TEST_BINARY_DIR}"
    "${TEST_DIR}/test.c"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${_KLU_INCLUDE_DIRS_ARG}"
    LINK_LIBRARIES ${_KLU_LINK_LIBRARIES}
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile KLU test. Check output in ${TEST_DIR}/compile.out")
  endif()

else()
  message(STATUS "Skipped KLU checks.")
endif()
