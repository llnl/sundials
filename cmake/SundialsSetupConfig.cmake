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
# Configures the SUNDIALS config header file:
#  sundials_config.h
# ---------------------------------------------------------------

# ============================================================================
# Generate macros and substitution variables related to the build options.
# ============================================================================

# prepare substitution variable PRECISION_LEVEL for sundials_config.h
string(TOUPPER ${SUNDIALS_PRECISION} SUNDIALS_PRECISION)
set(PRECISION_LEVEL "#define SUNDIALS_${SUNDIALS_PRECISION}_PRECISION 1")

# prepare substitution variable INDEX_TYPE for sundials_config.h
set(INDEX_TYPE "#define SUNDIALS_INT${SUNDIALS_INDEX_SIZE}_T 1")

if(COMPILER_HAS_DEPRECATED_MSG)
  set(SUNDIALS_DEPRECATED_MSG_MACRO "${COMPILER_DEPRECATED_MSG_ATTRIBUTE}")
else()
  set(SUNDIALS_DEPRECATED_MSG_MACRO "SUNDIALS_DEPRECATED")
endif()

if($ENV{CI_JOB_ID})
  set(JOB_ID $ENV{CI_JOB_ID})
else()
  string(TIMESTAMP JOB_ID "%Y%m%d%H%M%S")
endif()

if($ENV{CI_JOB_STARTED_AT})
  set(JOB_START_TIME $ENV{CI_JOB_STARTED_AT})
else()
  string(TIMESTAMP JOB_START_TIME "%Y%m%d%H%M%S")
endif()

# ============================================================================
# Generate macros and substitution variables related to TPLs that SUNDIALS is
# being built with.
# ============================================================================

# prepare substitution variables for modules that have been built
set(SUNDIALS_CONFIGH_BUILDS "")
foreach(_item ${SUNDIALS_BUILD_LIST})
  if(${${_item}})
    string(REPLACE "SUNDIALS_ENABLE_" "" _module ${_item})
    string(APPEND SUNDIALS_CONFIGH_BUILDS "#define SUNDIALS_${_module} 1\n")
  endif()
endforeach()

# prepare substitution variable SUNDIALS_${TPL NAME}_ENABLED for
# sundials_config.h
foreach(tpl ${SUNDIALS_TPL_LIST})
  set(SUNDIALS_${tpl}_ENABLED TRUE)
endforeach()

# prepare substitution variable SUNDIALS_TRILINOS_HAVE_MPI for sundials_config.h
if(SUNDIALS_ENABLE_MPI)
  set(SUNDIALS_TRILINOS_HAVE_MPI TRUE)
endif()

# prepare substitution variable(s) SUNDIALS_RAJA_BACKENDS_*
foreach(backend ${SUNDIALS_RAJA_BACKENDS})
  set(SUNDIALS_RAJA_BACKENDS_${backend} TRUE)
endforeach()

# prepare substitution variable(s) SUNDIALS_GINKGO_BACKENDS_*
foreach(backend ${SUNDIALS_GINKGO_BACKENDS})
  set(SUNDIALS_GINKGO_BACKENDS_${backend} TRUE)
endforeach()

# prepare substitution variable(s) SUNDIALS_MAGMA_BACKENDS_*
foreach(backend ${SUNDIALS_MAGMA_BACKENDS})
  set(SUNDIALS_MAGMA_BACKENDS_${backend} TRUE)
endforeach()

# prepare substitution variable SUNDIALS_HAVE_POSIX_TIMERS for sundials_config.h
if(SUNDIALS_POSIX_TIMERS) # set in SundialsPOSIXTimers.cmake
  set(SUNDIALS_HAVE_POSIX_TIMERS TRUE)
endif()

# =============================================================================
# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.
# =============================================================================

configure_file(${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
               ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h)
