# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# Module that sets up compilers for SUNDIALS.
# ---------------------------------------------------------------

# ===============================================================
# Determine the index type for the compiler
# ===============================================================

include(SundialsIndexSize)

# ===============================================================
# Platform specific settings
# ===============================================================

if(WIN32)
  # Under Windows, add compiler directive to inhibit warnings about use of
  # insecure functions.
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)

  # Under Windows, we need to have dll and exe files in the same directory to
  # run the test suite properly.
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

if(APPLE)
  # Allow undefined symbols that will be resolved by a user program.
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
endif()

# ===============================================================
# RPath settings
# ===============================================================

# only apply rpath settings for builds using shared libs
if(BUILD_SHARED_LIBS)
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already (but later on when
  # installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")

  # add the automatically determined parts of the RPATH which point to
  # directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system
  # directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
       "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  endif()
endif()

# ===============================================================
# Configure compiler flags
#
# TODO(DJG): Set flags based on CMAKE_<language>_COMPILER_ID
# ===============================================================

if(ENABLE_ALL_WARNINGS)
  message(STATUS "Enabling all compiler warnings")

  # Some warning flags are not supported by all compilers so ignore unknown
  # flags with -Wno-unknown-warning-option. Ironically, this is not supported by
  # some compilers.
  set(WARNING_FLAGS
      "-Wno-unknown-warning-option -Wall -Wpedantic -Wextra -Wshadow \
-Wwrite-strings -Wcast-align -Wdisabled-optimization -Wvla -Walloca \
-Wduplicated-cond -Wduplicated-branches -Wunused-macros \
-Wunused-local-typedefs")
  # TODO(SBR): Try to add -Wredundant-decls once SuperLU version is updated in
  # CI tests
  # TODO(DJG): To add -Wundef need to workaround __cplusplus and GTest warnings

  # Avoid numerous warnings from printf
  if(SUNDIALS_PRECISION MATCHES "EXTENDED")
    set(WARNING_FLAGS "-Wdouble-promotion ${WARNING_FLAGS}")
  endif()

  if((SUNDIALS_PRECISION MATCHES "DOUBLE") AND (SUNDIALS_INDEX_SIZE MATCHES "32"
                                               ))
    set(WARNING_FLAGS "-Wconversion -Wno-sign-conversion ${WARNING_FLAGS}")
  endif()

  # Avoid numerous warnings from SWIG generated functions
  if(NOT BUILD_FORTRAN_MODULE_INTERFACE)
    set(WARNING_FLAGS "-Wmissing-declarations -Wcast-qual ${WARNING_FLAGS}")
  endif()

  set(CMAKE_C_FLAGS "${WARNING_FLAGS} ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${WARNING_FLAGS} ${CMAKE_CXX_FLAGS}")

  # TODO(DJG): Add -fcheck=all,no-pointer,no-recursion once Jenkins is updated
  # to use gfortran > 5.5 which segfaults with -fcheck=array-temps,bounds,do,mem
  # no- options were added in gfortran 6
  #
  # Exclude run-time pointer checks (no-pointer) because passing null objects to
  # SUNDIALS functions (e.g., sunmat => null() to SetLinearSolver) causes a
  # run-time error with this check
  #
  # Exclude checks for subroutines and functions not marked as recursive
  # (no-recursion) e.g., ark_brusselator1D_task_local_nls_f2003 calls
  # SUNNonlinsolFree from within a custom nonlinear solver implementation of
  # SUNNonlinsolFree which causes a run-time error with this check
  set(CMAKE_Fortran_FLAGS
      "-Wall -Wpedantic -Wno-unused-dummy-argument -Wno-c-binding-type -ffpe-summary=none ${CMAKE_Fortran_FLAGS}"
  )
endif()

if(ENABLE_WARNINGS_AS_ERRORS)
  message(STATUS "Enabling compiler warnings as errors")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --Werror all-warnings") # CUDA 10.2+
endif()

# With clang it is not possible to combine the -fsanitize=address and
# -fsanitize=memory checkers.

if(ENABLE_ADDRESS_SANITIZER)
  message(STATUS "Enabling address sanitizer")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=address")
endif()

if(ENABLE_MEMORY_SANITIZER)
  message(STATUS "Enabling memory sanitizer")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=memory")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=memory")
endif()

if(ENABLE_LEAK_SANITIZER)
  message(STATUS "Enabling leak sanitizer")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=leak")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=leak")
endif()

if(ENABLE_UNDEFINED_BEHAVIOR_SANITIZER)
  message(STATUS "Enabling undefined behavior sanitizer")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=undefined")
endif()

if(SUNDIALS_DEBUG)
  message(STATUS "Adding debugging preprocessor directives")

  foreach(debug ${_SUNDIALS_DEBUG_OPTIONS})
    if(${${debug}})
      add_compile_definitions(${debug})
    endif()
  endforeach()
endif()

# ===============================================================
# C settings
# ===============================================================

set(DOCSTR "The C standard to use (99, 11, 17, 23)")
sundials_option(CMAKE_C_STANDARD STRING "${DOCSTR}" "99" OPTIONS "99;11;17;23")
message(STATUS "C standard set to ${CMAKE_C_STANDARD}")

set(DOCSTR "Enable C compiler specific extensions")
sundials_option(CMAKE_C_EXTENSIONS BOOL "${DOCSTR}" ON)
message(STATUS "C extensions set to ${CMAKE_C_EXTENSIONS}")

# ---------------------------------------------------------------
# Check for __builtin_expect
# ---------------------------------------------------------------

check_c_source_compiles(
  "
  #include <stdio.h>
  int main(void) {
    double a = 0.0;
    if (__builtin_expect(a < 0, 0)) {
      a = 0.0;
    }
    a = a + 1.0;
    printf(\"a=%g\", a);
    return 0;
  }
"
  SUNDIALS_C_COMPILER_HAS_BUILTIN_EXPECT)

# ---------------------------------------------------------------
# Check for assume related extensions
# ---------------------------------------------------------------

# gcc >= 13 should have __attribute__((assume))
check_c_source_compiles(
  "
  #include <stdio.h>
  int main(void) {
    double a = 0.0;
    #if defined(__has_attribute)
    # if !__has_attribute(assume)
    #   error no assume
    # endif
    #else
    #error no __has_attribute
    #endif
    __attribute__((assume(a >= 0.0)));
    a = a + 1.0;
    printf(\"a=%g\", a);
    return 0;
  }
"
  SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_ASSUME)

# LLVM based compilers should have __builtin_assume
if(NOT SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_ASSUME)
  check_c_source_compiles(
    "
    #include <stdio.h>
    int main(void) {
      double a = 0.0;
      __builtin_assume(a >= 0.0);
      a = a + 1.0;
      printf(\"a=%g\", a);
      return 0;
    }
  "
    SUNDIALS_C_COMPILER_HAS_BUILTIN_ASSUME)
endif()

# MSVC provides __assume
if(NOT (SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_ASSUME
        OR SUNDIALS_C_COMPILER_HAS_BUILTIN_ASSUME))
  check_c_source_compiles(
    "
    #include <stdio.h>
    int main(void) {
      double a = 0.0;
      __assume(a >= 0.0));
      a = a + 1.0;
      printf(\"a=%g\", a);
      return 0;
    }
  "
    SUNDIALS_C_COMPILER_HAS_ASSUME)
endif()

# ---------------------------------------------------------------
# Check for unused extension
# ---------------------------------------------------------------

check_c_source_compiles(
  "
  int main(void) {
    __attribute__((unused)) double a = 0.0;
    return 0;
  }
"
  SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_UNUSED)

# ---------------------------------------------------------------
# Check for POSIX timers
# ---------------------------------------------------------------
include(SundialsPOSIXTimers)

if(SUNDIALS_POSIX_TIMERS AND POSIX_TIMERS_NEED_POSIX_C_SOURCE)
  set(DOCSTR "Value of _POSIX_C_SOURCE")
  sundials_option(SUNDIALS_POSIX_C_SOURCE STRING "${DOCSTR}" "200112L" ADVANCED)
  set(CMAKE_C_FLAGS
      "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=${SUNDIALS_POSIX_C_SOURCE}")
endif()

# ---------------------------------------------------------------
# Check for deprecated attribute with message
# ---------------------------------------------------------------
if(WIN32)
  set(COMPILER_DEPRECATED_MSG_ATTRIBUTE
      "__declspec(deprecated(msg))"
      CACHE INTERNAL "")
else()
  set(COMPILER_DEPRECATED_MSG_ATTRIBUTE
      "__attribute__ ((__deprecated__(msg)))"
      CACHE INTERNAL "")
endif()
check_c_source_compiles(
  "
  #define msg \"test\"
  ${COMPILER_DEPRECATED_MSG_ATTRIBUTE} int somefunc(void) { return 0; }
  int main(void) { return somefunc();}"
  COMPILER_HAS_DEPRECATED_MSG)

# ===============================================================
# Fortran settings
# ===============================================================

# When LAPACK is enabled we will need a Fortran compiler to infer the
# name-mangling scheme if it is not set by the user
if(ENABLE_LAPACK)
  set(NEED_FORTRAN_NAME_MANGLING TRUE)
endif()

# ------------------------------------------------------------------------------
# Allow the user to manually specify the Fortran name-mangling scheme
#
# The build system tries to infer the Fortran name-mangling scheme using a
# Fortran compiler and defaults to using lower case and one underscore if the
# scheme can not be determined. If a working Fortran compiler is not available
# or the user needs to override the inferred or default scheme, the following
# options specify the case and number of appended underscores corresponding to
# the Fortran name-mangling scheme of symbol names that do not themselves
# contain underscores. This is all we really need for the LAPACK interfaces. A
# working Fortran compiler is only necessary for building Fortran example
# programs.
# ------------------------------------------------------------------------------

# The case to use in the name-mangling scheme
sundials_option(SUNDIALS_LAPACK_CASE STRING
                "case of LAPACK function names (lower/upper)" "" ADVANCED)

# The number of underscores of appended in the name-mangling scheme
sundials_option(
  SUNDIALS_LAPACK_UNDERSCORES STRING
  "number of underscores appended to LAPACK function names (none/one/two)" ""
  ADVANCED)

# If used, both case and underscores must be set
if((NOT SUNDIALS_LAPACK_CASE) AND SUNDIALS_LAPACK_UNDERSCORES)
  message(FATAL_ERROR "If SUNDIALS_LAPACK_UNDERSCORES is set, "
                      "SUNDIALS_LAPACK_CASE must also be set.")
endif()
if(SUNDIALS_LAPACK_CASE AND (NOT SUNDIALS_LAPACK_UNDERSCORES))
  message(FATAL_ERROR "If SUNDIALS_LAPACK_CASE is set, "
                      "SUNDIALS_LAPACK_UNDERSCORES must also be set.")
endif()

# Did the user provide a name-mangling scheme?
if(SUNDIALS_LAPACK_CASE AND SUNDIALS_LAPACK_UNDERSCORES)

  string(TOUPPER ${SUNDIALS_LAPACK_CASE} SUNDIALS_LAPACK_CASE)
  string(TOUPPER ${SUNDIALS_LAPACK_UNDERSCORES} SUNDIALS_LAPACK_UNDERSCORES)

  # Based on the given case and number of underscores, set the C preprocessor
  # macro definitions. Since SUNDIALS never uses symbols names containing
  # underscores we set the name-mangling schemes to be the same. In general,
  # names of symbols with and without underscore may be mangled differently
  # (e.g. g77 mangles mysub to mysub_ and my_sub to my_sub__)
  if(SUNDIALS_LAPACK_CASE MATCHES "LOWER")
    if(SUNDIALS_LAPACK_UNDERSCORES MATCHES "NONE")
      set(LAPACK_MANGLE_MACRO1 "#define SUNDIALS_LAPACK_FUNC(name,NAME) name")
      set(LAPACK_MANGLE_MACRO2 "#define SUNDIALS_LAPACK_FUNC_(name,NAME) name")
    elseif(SUNDIALS_LAPACK_UNDERSCORES MATCHES "ONE")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) name ## _")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) name ## _")
    elseif(SUNDIALS_LAPACK_UNDERSCORES MATCHES "TWO")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) name ## __")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) name ## __")
    else()
      message(FATAL_ERROR "Invalid SUNDIALS_LAPACK_UNDERSCORES option.")
    endif()
  elseif(SUNDIALS_LAPACK_CASE MATCHES "UPPER")
    if(SUNDIALS_LAPACK_UNDERSCORES MATCHES "NONE")
      set(LAPACK_MANGLE_MACRO1 "#define SUNDIALS_LAPACK_FUNC(name,NAME) NAME")
      set(LAPACK_MANGLE_MACRO2 "#define SUNDIALS_LAPACK_FUNC_(name,NAME) NAME")
    elseif(SUNDIALS_LAPACK_UNDERSCORES MATCHES "ONE")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) NAME ## _")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) NAME ## _")
    elseif(SUNDIALS_LAPACK_UNDERSCORES MATCHES "TWO")
      set(LAPACK_MANGLE_MACRO1
          "#define SUNDIALS_LAPACK_FUNC(name,NAME) NAME ## __")
      set(LAPACK_MANGLE_MACRO2
          "#define SUNDIALS_LAPACK_FUNC_(name,NAME) NAME ## __")
    else()
      message(FATAL_ERROR "Invalid SUNDIALS_LAPACK_UNDERSCORES option.")
    endif()
  else()
    message(FATAL_ERROR "Invalid SUNDIALS_LAPACK_CASE option.")
  endif()

  # name-mangling scheme has been manually set
  set(NEED_FORTRAN_NAME_MANGLING FALSE)

  configure_file(${PROJECT_SOURCE_DIR}/src/sundials/sundials_lapack_defs.h.in
                 ${PROJECT_BINARY_DIR}/src/sundials/sundials_lapack_defs.h)

endif()

# Do we need a Fortran compiler?
if(BUILD_FORTRAN_MODULE_INTERFACE OR NEED_FORTRAN_NAME_MANGLING)
  include(SundialsSetupFortran)
endif()

# ===============================================================
# C++ settings
# ===============================================================

if(BUILD_BENCHMARKS
   OR SUNDIALS_TEST_ENABLE_UNIT_TESTS
   OR EXAMPLES_ENABLE_CXX
   OR ENABLE_CUDA
   OR ENABLE_HIP
   OR ENABLE_SYCL
   OR ENABLE_RAJA
   OR ENABLE_TRILINOS
   OR ENABLE_SUPERLUDIST
   OR ENABLE_MAGMA
   OR ENABLE_GINKGO
   OR ENABLE_KOKKOS
   OR ENABLE_ADIAK)
  include(SundialsSetupCXX)
endif()

# ===============================================================
# CUDA settings
# ===============================================================

if(ENABLE_CUDA)
  include(SundialsSetupCuda)
  # we treat CUDA as both a TPL and a language
  list(APPEND SUNDIALS_TPL_LIST "CUDA")
endif()

# ===============================================================
# HIP settings
# ===============================================================

if(ENABLE_HIP)
  include(SundialsSetupHIP)
  # we treat HIP as both a TPL and a language
  list(APPEND SUNDIALS_TPL_LIST "HIP")
endif()

# ===============================================================
# Configure presentation of language options
# ===============================================================

# List of enabled languages
set(_SUNDIALS_ENABLED_LANGS "C")
if(CXX_FOUND)
  list(APPEND _SUNDIALS_ENABLED_LANGS "CXX")
endif()
if(Fortran_FOUND)
  list(APPEND _SUNDIALS_ENABLED_LANGS "Fortran")
endif()
if(CUDA_FOUND)
  list(APPEND _SUNDIALS_ENABLED_LANGS "CUDA")
endif()

# Upper case version of build type
string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)

# Make build type specific flag options ADVANCED, except for the one
# corresponding to the current build type
foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  foreach(build_type DEBUG;RELEASE;RELWITHDEBINFO;MINSIZEREL)
    if("${_cmake_build_type}" STREQUAL "${build_type}")
      message(STATUS "Appending ${lang} ${build_type} flags")
      mark_as_advanced(CLEAR CMAKE_${lang}_FLAGS_${build_type})
    else()
      mark_as_advanced(FORCE CMAKE_${lang}_FLAGS_${build_type})
    endif()
  endforeach()
  # show the language compiler and flags
  mark_as_advanced(CLEAR CMAKE_${lang}_COMPILER CMAKE_${lang}_FLAGS)
endforeach()

# ===============================================================
# Configure compilers for installed examples
# ===============================================================

foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  if(ENABLE_MPI)
    if(DEFINED MPI_${lang}_COMPILER)
      set(_EXAMPLES_${lang}_COMPILER
          "${MPI_${lang}_COMPILER}"
          CACHE INTERNAL "${lang} compiler for installed examples")
    endif()
  else()
    set(_EXAMPLES_${lang}_COMPILER
        "${CMAKE_${lang}_COMPILER}"
        CACHE INTERNAL "${lang} compiler for installed examples")
  endif()
endforeach()

# ===============================================================
# Configure clang-tidy for linting
# ===============================================================

set(SUNDIALS_DEV_CLANG_TIDY_DIR ${CMAKE_BINARY_DIR}/clang-tidy/)

if(SUNDIALS_DEV_CLANG_TIDY)
  find_program(CLANG_TIDY_PATH NAMES clang-tidy)
  if(NOT CLANG_TIDY_PATH)
    message(FATAL_ERROR "Could not find the program clang-tidy")
  endif()
  message(STATUS "Found clang-tidy: ${CLANG_TIDY_PATH}")

  make_directory(${SUNDIALS_DEV_CLANG_TIDY_DIR})
  if(SUNDIALS_DEV_CLANG_TIDY_FIX_ERRORS)
    set(CMAKE_C_CLANG_TIDY ${CLANG_TIDY_PATH} -format-style='file' --fix)
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_PATH} -format-style='file' --fix)
  else()
    set(CMAKE_C_CLANG_TIDY
        ${CLANG_TIDY_PATH} -format-style='file'
        --export-fixes=${SUNDIALS_DEV_CLANG_TIDY_DIR}/clang-tidy-fixes.yaml)
    set(CMAKE_CXX_CLANG_TIDY
        ${CLANG_TIDY_PATH} -format-style='file'
        --export-fixes=${SUNDIALS_DEV_CLANG_TIDY_DIR}/clang-tidy-cxx-fixes.yaml)
  endif()
endif()

if(SUNDIALS_DEV_IWYU)
  find_program(IWYU_PATH NAMES include-what-you-use iwyu)
  if(NOT IWYU_PATH)
    message(FATAL_ERROR "Could not find the program include-what-you-use")
  endif()
  message(STATUS "Found IWYU: ${IWYU_PATH}")
  set(CMAKE_C_INCLUDE_WHAT_YOU_USE
      ${IWYU_PATH} -Xiwyu --mapping_file=${CMAKE_SOURCE_DIR}/scripts/iwyu.imp
      -Xiwyu --error_always)
  set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE
      ${IWYU_PATH} -Xiwyu --mapping_file=${CMAKE_SOURCE_DIR}/scripts/iwyu.imp
      -Xiwyu --error_always)
endif()
