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

#[=======================================================================[.rst:
SundialsAddWarningFlags
-----------------------

This module provides a command for adding compiler warning flags.

Load this module with:

.. code-block:: cmake

   include(SundialsAddWarningFlags)

Commands
^^^^^^^^

This module provides the following command:

.. cmake:command:: sundials_add_warning_flags

   Add compiler warning flags for a language.

   .. code-block:: cmake

      sundials_add_warning_flags(<lang>)

   Prepends warning flags to ``CMAKE_<lang>_FLAGS`` so users can override them.
   The flags are selected based on the compiler ID and command-line style
   (GCC-like or MSVC-like), so call this command after ``<lang>`` is enabled.

   The arguments are:

   ``<lang>``
     The language to add flags for: ``C``, ``CXX``, ``CUDA``, or ``Fortran``.

   The flags added depend on the following options:

   :cmakeop:`SUNDIALS_ENABLE_ALL_WARNINGS`
     Add SUNDIALS warning flags for C, C++, and GNU Fortran compilers.

     * Sets ``SUNDIALS_<lang>_SWIG_WARNING_FLAGS`` to a list of flags that
       disable warnings in SWIG generated code that cannot be fixed. These
       flags are used by ``sundials_add_f2003_library``.

   :cmakeop:`CMAKE_COMPILE_WARNING_AS_ERROR`
     Add flags to treat warnings as errors with CMake versions prior to 3.24.

     * Newer versions of CMake add these flags to targets natively.

     * This option is independent of :cmakeop:`SUNDIALS_ENABLE_ALL_WARNINGS`.
#]=======================================================================]

function(sundials_add_warning_flags lang)
  set(_id "${CMAKE_${lang}_COMPILER_ID}")

  set(_msvc_style FALSE)
  if(_id STREQUAL "MSVC" OR CMAKE_${lang}_COMPILER_FRONTEND_VARIANT STREQUAL
                            "MSVC")
    set(_msvc_style TRUE)
  endif()

  set(_clang_ids
      Clang
      AppleClang
      ARMClang
      CrayClang
      FujitsuClang
      IBMClang
      IntelLLVM)

  set(_flags)
  set(_swig_flags)

  if(SUNDIALS_ENABLE_ALL_WARNINGS)
    if(lang STREQUAL "C" OR lang STREQUAL "CXX")
      if(_msvc_style)
        # MSVC-style command line (MSVC, clang-cl, icx-cl)
        list(APPEND _flags /W4)
      else()
        # GCC-style command line
        list(
          APPEND
          _flags
          -Wall
          -Wpedantic
          -Wextra
          -Wshadow
          -Wwrite-strings
          -Wcast-align
          -Wcast-qual
          -Wdisabled-optimization
          -Wvla
          -Walloca
          -Wduplicated-cond
          -Wduplicated-branches
          -Wmissing-declarations
          -Wunused-macros
          -Wunused-local-typedefs
          -Wundef)
        # TODO(SBR): Try to add -Wredundant-decls once SuperLU version is
        # updated in CI tests

        # Avoid numerous warnings from printf
        if(SUNDIALS_PRECISION MATCHES "EXTENDED")
          list(APPEND _flags -Wdouble-promotion)
        endif()

        # -Wno-sign-conversion must follow -Wconversion (Clang applies the flags
        # in order, so the reverse re-enables sign conversion warnings)
        if(SUNDIALS_PRECISION MATCHES "DOUBLE" AND SUNDIALS_INDEX_SIZE MATCHES
                                                   "32")
          list(APPEND _flags -Wconversion -Wno-sign-conversion)
        endif()

        list(APPEND _swig_flags -Wno-cast-qual -Wno-missing-declarations
             -Wno-unused-macros)
      endif()

      # Clang-based compilers: ignore unsupported flags (GCC-only or newer than
      # the compiler) and add Clang-only warnings. GCC rejects unknown -W flags,
      # so these are not added for other compilers.
      if(_id IN_LIST _clang_ids)
        list(APPEND _flags -Wno-unknown-warning-option -Wreserved-identifier)
        list(APPEND _swig_flags -Wno-reserved-identifier)
      endif()
    elseif(lang STREQUAL "Fortran" AND _id STREQUAL "GNU")
      # TODO(DJG): Add -fcheck=all,no-pointer,no-recursion once Jenkins is
      # updated to use gfortran > 5.5 which segfaults with
      # -fcheck=array-temps,bounds,do,mem no- options were added in gfortran 6
      #
      # Exclude run-time pointer checks (no-pointer) because passing null
      # objects to SUNDIALS functions (e.g., sunmat => null() to
      # SetLinearSolver) causes a run-time error with this check
      #
      # Exclude checks for subroutines and functions not marked as recursive
      # (no-recursion) e.g., ark_brusselator1D_task_local_nls_f2003 calls
      # SUNNonlinsolFree from within a custom nonlinear solver implementation of
      # SUNNonlinsolFree which causes a run-time error with this check
      list(
        APPEND
        _flags
        -Wall
        -Wpedantic
        -Wno-unused-dummy-argument
        -Wno-c-binding-type
        -ffpe-summary=none)
    endif()
  endif()

  if(CMAKE_COMPILE_WARNING_AS_ERROR AND CMAKE_VERSION VERSION_LESS 3.24)
    if(lang STREQUAL "Fortran" AND _id MATCHES "^Intel")
      if(_msvc_style)
        list(APPEND _flags -warn:errors)
      else()
        list(APPEND _flags -warn errors)
      endif()
    elseif(lang STREQUAL "CUDA" AND _id STREQUAL "NVIDIA")
      list(APPEND _flags -Werror all-warnings) # CUDA 10.2+
    elseif(_msvc_style)
      list(APPEND _flags -WX)
    else()
      list(APPEND _flags -Werror)
    endif()
  endif()

  # Prepend our flags so users can override them with CMAKE_<lang>_FLAGS
  if(_flags)
    list(JOIN _flags " " _flags)
    string(STRIP "${_flags} ${CMAKE_${lang}_FLAGS}" _flags)
    set(CMAKE_${lang}_FLAGS
        "${_flags}"
        PARENT_SCOPE)
  endif()
  set(SUNDIALS_${lang}_SWIG_WARNING_FLAGS
      "${_swig_flags}"
      PARENT_SCOPE)
endfunction()
