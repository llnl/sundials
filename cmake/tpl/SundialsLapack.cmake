# -----------------------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
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
# Module to find and setup LAPACK/BLAS correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# <TPL>_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# LAPACK does not support extended precision
if(ENABLE_LAPACK AND SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR "LAPACK is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(LAPACK REQUIRED)

# get path to LAPACK library to use in generated makefiles for examples, if
# LAPACK_LIBRARIES contains multiple items only use the path of the first entry
list(GET LAPACK_LIBRARIES 0 TMP_LAPACK_LIBRARIES)
get_filename_component(LAPACK_LIBRARY_DIR ${TMP_LAPACK_LIBRARIES} PATH)

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# ---------------------------------------------------------------
# Determining the name-mangling scheme if needed
# ---------------------------------------------------------------
# In general, names of symbols with and without underscore may be mangled
# differently (e.g. g77 mangles mysub to mysub_ and my_sub to my_sub__), we have
# to consider both cases.
#
# Method:
#
# 1. create a library from a Fortran source file which defines a function "mysub"
# 2. attempt to link with this library a C source file which calls the "mysub"
#    function using various possible schemes (6 different schemes, corresponding
#    to all combinations lower/upper case and none/one/two underscores).
# 3. define the name-mangling scheme based on the test that was successful.
#
# On exit, if we were able to infer the scheme, the variables
# CMAKE_Fortran_SCHEME_NO_UNDERSCORES and CMAKE_Fortran_SCHEME_WITH_UNDERSCORES
# contain the mangled names for "mysub" and "my_sub", respectively.
# ---------------------------------------------------------------
if(NEED_FORTRAN_NAME_MANGLING)

  set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "")
  set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES "")

  # Create the FortranTest directory
  set(FortranTest_DIR ${PROJECT_BINARY_DIR}/FortranTest)
  file(MAKE_DIRECTORY ${FortranTest_DIR})

  # Create a CMakeLists.txt file which will generate the "flib" library and an
  # executable "ftest"
  file(
    WRITE ${FortranTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ftest Fortran)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
    "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
    "ADD_LIBRARY(flib flib.f)\n"
    "ADD_EXECUTABLE(ftest ftest.f)\n"
    "TARGET_LINK_LIBRARIES(ftest flib)\n")

  # Create the Fortran source flib.f which defines two subroutines, "mysub" and
  # "my_sub"
  file(WRITE ${FortranTest_DIR}/flib.f
       "        SUBROUTINE mysub\n" "        RETURN\n" "        END\n"
       "        SUBROUTINE my_sub\n" "        RETURN\n" "        END\n")

  # Create the Fortran source ftest.f which calls "mysub" and "my_sub"
  file(WRITE ${FortranTest_DIR}/ftest.f
       "        PROGRAM ftest\n" "        CALL mysub()\n"
       "        CALL my_sub()\n" "        END\n")

  # Use TRY_COMPILE to make the targets "flib" and "ftest"
  try_compile(
    FTEST_OK ${FortranTest_DIR}
    ${FortranTest_DIR} ftest
    OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, we must remove the
  # CMakeFiles directory.
  file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)

  # Proceed based on test results
  if(FTEST_OK)

    # Infer Fortran name-mangling scheme for symbols WITHOUT underscores.
    # Overwrite CMakeLists.txt with one which will generate the "ctest1"
    # executable
    file(
      WRITE ${FortranTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
      "PROJECT(ctest1 C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ctest1 ctest1.c)\n"
      "FIND_LIBRARY(FLIB flib \"${FortranTest_DIR}\")\n"
      "TARGET_LINK_LIBRARIES(ctest1 \${FLIB})\n")

    # Define the list "options" of all possible schemes that we want to consider
    # Get its length and initialize the counter "iopt" to zero
    set(options mysub mysub_ mysub__ MYSUB MYSUB_ MYSUB__)
    list(LENGTH options imax)
    set(iopt 0)

    # We will attempt to successfully generate the "ctest1" executable as long
    # as there still are entries in the "options" list
    while(${iopt} LESS ${imax})
      # Get the current list entry (current scheme)
      list(GET options ${iopt} opt)
      # Generate C source which calls the "mysub" function using the current
      # scheme
      file(WRITE ${FortranTest_DIR}/ctest1.c
           "extern void ${opt}();\n" "int main(void){${opt}();return(0);}\n")
      # Use TRY_COMPILE to make the "ctest1" executable from the current C
      # source and linking to the previously created "flib" library.
      try_compile(
        CTEST_OK ${FortranTest_DIR}
        ${FortranTest_DIR} ctest1
        OUTPUT_VARIABLE MY_OUTPUT)
      # Write output compiling the test code
      file(WRITE ${FortranTest_DIR}/ctest1_${opt}.out "${MY_OUTPUT}")
      # To ensure we do not use stuff from the previous attempts, we must remove
      # the CMakeFiles directory.
      file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
      # Test if we successfully created the "ctest" executable. If yes, save the
      # current scheme, and set the counter "iopt" to "imax" so that we exit the
      # while loop. Otherwise, increment the counter "iopt" and go back in the
      # while loop.
      if(CTEST_OK)
        set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES ${opt})
        set(iopt ${imax})
      else(CTEST_OK)
        math(EXPR iopt ${iopt}+1)
      endif()
    endwhile(${iopt} LESS ${imax})

    # Infer Fortran name-mangling scheme for symbols WITH underscores.
    # Practically a duplicate of the previous steps.
    file(
      WRITE ${FortranTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
      "PROJECT(ctest2 C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ctest2 ctest2.c)\n"
      "FIND_LIBRARY(FLIB flib \"${FortranTest_DIR}\")\n"
      "TARGET_LINK_LIBRARIES(ctest2 \${FLIB})\n")

    set(options my_sub my_sub_ my_sub__ MY_SUB MY_SUB_ MY_SUB__)
    list(LENGTH options imax)
    set(iopt 0)
    while(${iopt} LESS ${imax})
      list(GET options ${iopt} opt)
      file(WRITE ${FortranTest_DIR}/ctest2.c
           "extern void ${opt}();\n" "int main(void){${opt}();return(0);}\n")
      try_compile(
        CTEST_OK ${FortranTest_DIR}
        ${FortranTest_DIR} ctest2
        OUTPUT_VARIABLE MY_OUTPUT)
      file(WRITE ${FortranTest_DIR}/ctest2_${opt}.out "${MY_OUTPUT}")
      file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
      if(CTEST_OK)
        set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES ${opt})
        set(iopt ${imax})
      else(CTEST_OK)
        math(EXPR iopt ${iopt}+1)
      endif()
    endwhile(${iopt} LESS ${imax})

    # If a name-mangling scheme was found set the C preprocessor macros to use
    # that scheme. Otherwise default to lower case with one underscore.
    if(CMAKE_Fortran_SCHEME_NO_UNDERSCORES
       AND CMAKE_Fortran_SCHEME_WITH_UNDERSCORES)
      message(STATUS "Determining Fortran name-mangling scheme... OK")
    else()
      message(STATUS "Determining Fortran name-mangling scheme... DEFAULT")
      set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "mysub_")
      set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES "my_sub_")
    endif()

    # Symbols NO underscores
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub")
      set(LAPACK_MANGLE_MACRO1 "#define SUNDIALS_LAPACK_FUNC(name,NAME) name")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub_")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) name ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub__")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) name ## __")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB")
      set(LAPACK_MANGLE_MACRO1 "#define SUNDIALS_LAPACK_FUNC(name,NAME) NAME")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB_")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) NAME ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB__")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) NAME ## __")
    endif()

    # Symbols WITH underscores
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "my_sub")
      set(LAPACK_MANGLE_MACRO2 "#define SUNDIALS_LAPACK_FUNC_(name,NAME) name")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "my_sub_")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) name ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "my_sub__")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) name ## __")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "MY_SUB")
      set(LAPACK_MANGLE_MACRO2 "#define SUNDIALS_LAPACK_FUNC_(name,NAME) NAME")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "MY_SUB_")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) NAME ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "MY_SUB__")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) NAME ## __")
    endif()

    # name-mangling scheme has been set
    set(NEED_FORTRAN_NAME_MANGLING FALSE)

    configure_file(${PROJECT_SOURCE_DIR}/src/sundials/sundials_lapack_defs.h.in
                   ${PROJECT_BINARY_DIR}/src/sundials/sundials_lapack_defs.h)

  else(FTEST_OK)
    message(STATUS "Determining Fortran name-mangling scheme... FAILED")
  endif()

endif()

# Try building a simple test
if(NOT LAPACK_WORKS)

  message(CHECK_START "Testing LAPACK")

  # Create the test directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/LAPACK_TEST)

  # Create a C source file calling a BLAS (dcopy) and LAPACK (dgetrf) function
  file(
    WRITE ${TEST_DIR}/test.c
    "${LAPACK_MANGLE_MACRO1}\n"
    "#define dcopy_f77 SUNDIALS_LAPACK_FUNC(dcopy, DCOPY)\n"
    "#define dgetrf_f77 SUNDIALS_LAPACK_FUNC(dgetrf, DGETRF)\n"
    "extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);\n"
    "extern void dgetrf_f77(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);\n"
    "int main(void) {\n"
    "int n=1;\n"
    "double x=1.0;\n"
    "double y=1.0;\n"
    "dcopy_f77(&n, &x, &n, &y, &n);\n"
    "dgetrf_f77(&n, &n, &x, &n, &n, &n);\n"
    "return 0;\n"
    "}\n")

  # Workaround bug in older versions of CMake where the BLAS::BLAS target, which
  # LAPACK::LAPACK depends on, is not defined in the file
  # ${TEST_DIR}/CMakeFiles/CMakeTmp/<random_name>Targets.cmake created by
  # try_compile
  set(_lapack_targets LAPACK::LAPACK)
  if(CMAKE_VERSION VERSION_LESS 3.20)
    list(APPEND _lapack_targets BLAS::BLAS)
  endif()

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TEST_DIR}
    ${TEST_DIR}/test.c
    LINK_LIBRARIES ${_lapack_targets}
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Check the result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile LAPACK test. Check output in ${TEST_DIR}/compile.out"
    )
  endif()

else()
  message(STATUS "Skipped LAPACK test. Set LAPACK_WORKS=FALSE to test.")
endif()
