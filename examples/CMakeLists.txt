# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# We need this to ensure the installed templates have MPI when ENABLE_MPI=TRUE,
# at least until we convert all of the templates to use the SUNDIALS CMake
# target.

# ===================================================================
# Configure compilers for installed examples
# ===================================================================

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

# Disable some warnings for examples
if(ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif()

# Initialize CMAKE_Fortran_MODULE_DIRECTORY for examples
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ===================================================================
# Set variables used in generating CMake and Makefiles for examples
# ===================================================================

if(EXAMPLES_INSTALL)

  set(SHELL "sh")
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(exec_prefix "${CMAKE_INSTALL_PREFIX}")
  set(includedir "${prefix}/include")
  set(libdir "${exec_prefix}/${CMAKE_INSTALL_LIBDIR}")
  set(CPP "${CMAKE_CXX_COMPILER}")
  set(CPPFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CC "${CMAKE_C_COMPILER}")
  set(CFLAGS "${CMAKE_C_FLAGS_RELEASE}")
  set(LDFLAGS "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
  list2string(EXE_EXTRA_LINK_LIBS LIBS)
  list2string(EXE_EXTRA_LINK_LIBS EXAMPLES_EXTRA_LIBS)

  if(CXX_FOUND)
    set(CXX "${_EXAMPLES_CXX_COMPILER}")
    set(CXX_LNKR "${_EXAMPLES_CXX_COMPILER}")
    set(CXXFLAGS
        "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION}"
    )
    set(CXX_LDFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    list2string(EXE_EXTRA_LINK_LIBS CXX_LIBS)
  endif()

  if(SUPERLUMT_FOUND)
    list2string(SUPERLUMT_LIBRARIES SUPERLUMT_LIBS)
    set(SUPERLUMT_LIBS "${SUPERLUMT_LINKER_FLAGS} ${SUPERLUMT_LIBS}")
  endif()

  if(KLU_FOUND)
    list2string(KLU_LIBRARIES KLU_LIBS)
    set(KLU_LIBS "${KLU_LINKER_FLAGS} ${KLU_LIBS}")
  endif()

  if(LAPACK_FOUND)
    list2string(LAPACK_LIBRARIES LAPACK_LIBS)
  endif()

  if(HYPRE_FOUND)
    set(HYPRE_INC_DIR "${HYPRE_INCLUDE_DIR}")
    set(HYPRE_LIB_DIR "${HYPRE_LIBRARY_DIR}")
    set(HYPRE_LIBS "${HYPRE_LIBRARIES}")
  endif()

  if(MPI_C_FOUND)
    set(MPICC "${MPI_C_COMPILER}")
    set(MPI_INC_DIR ".")
    set(MPI_LIB_DIR ".")
    set(MPI_LIBS "")
    set(MPI_FLAGS "")
  endif()

  if(MPI_CXX_FOUND)
    set(MPICXX "${MPI_CXX_COMPILER}")
  endif()

endif()

# ===================================================================
# Add specific examples
# ===================================================================

# Add ARKode examples
if(BUILD_ARKODE)
  add_subdirectory(arkode)
endif()

# Add CVODE examples
if(BUILD_CVODE)
  add_subdirectory(cvode)
endif()

# Add CVODES Examples
if(BUILD_CVODES)
  add_subdirectory(cvodes)
endif()

# Add IDA examples
if(BUILD_IDA)
  add_subdirectory(ida)
endif()

# Add IDAS examples
if(BUILD_IDAS)
  add_subdirectory(idas)
endif()

# Add KINSOL examples
if(BUILD_KINSOL)
  add_subdirectory(kinsol)
endif()
