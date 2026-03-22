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

if(NOT
   (SUNDIALS_ENABLE_ARKODE
    AND SUNDIALS_ENABLE_CVODES
    AND SUNDIALS_ENABLE_IDAS
    AND SUNDIALS_ENABLE_KINSOL))
  message(
    FATAL_ERROR
      "sundials4py requires ARKODE, CVODES, IDAS, and KINSOL to be enabled")
endif()

# Warn if the user invokes CMake directly
if(NOT SKBUILD)
  message(
    WARNING
      "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

# nanobind needs the Python Interpreter and Development components
find_package(
  Python 3.12
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
  REQUIRED)

# Determine location of nanobind cmake config file
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)

# nanobind must already be installed as a Python module (e.g., with pip)
find_package(nanobind CONFIG REQUIRED)

# Add the source files for the bindings
set(sundials_SOURCES
    arkode/arkode_arkstep.cpp
    arkode/arkode_erkstep.cpp
    arkode/arkode_forcingstep.cpp
    arkode/arkode_lsrkstep.cpp
    arkode/arkode_mristep.cpp
    arkode/arkode_splittingstep.cpp
    arkode/arkode_sprkstep.cpp
    arkode/arkode.cpp
    cvodes/cvodes.cpp
    idas/idas.cpp
    kinsol/kinsol.cpp
    nvector/nvector_serial.cpp
    nvector/nvector_manyvector.cpp
    sundials4py.cpp
    sunadaptcontroller/sunadaptcontroller_imexgus.cpp
    sunadaptcontroller/sunadaptcontroller_mrihtol.cpp
    sunadaptcontroller/sunadaptcontroller_soderlind.cpp
    sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.cpp
    sundials/sundials_adaptcontroller.cpp
    sundials/sundials_adjointcheckpointscheme.cpp
    sundials/sundials_adjointstepper.cpp
    sundials/sundials_context.cpp
    sundials/sundials_core.cpp
    sundials/sundials_domeigestimator.cpp
    sundials/sundials_linearsolver.cpp
    sundials/sundials_logger.cpp
    sundials/sundials_matrix.cpp
    sundials/sundials_memory.cpp
    sundials/sundials_nonlinearsolver.cpp
    sundials/sundials_nvector.cpp
    sundials/sundials_profiler.cpp
    sundials/sundials_stepper.cpp
    sundomeigest/sundomeigest_power.cpp
    sunlinsol/sunlinsol_band.cpp
    sunlinsol/sunlinsol_dense.cpp
    sunlinsol/sunlinsol_pcg.cpp
    sunlinsol/sunlinsol_spbcgs.cpp
    sunlinsol/sunlinsol_spfgmr.cpp
    sunlinsol/sunlinsol_spgmr.cpp
    sunlinsol/sunlinsol_sptfqmr.cpp
    sunmatrix/sunmatrix_band.cpp
    sunmatrix/sunmatrix_dense.cpp
    sunmatrix/sunmatrix_sparse.cpp
    sunmemory/sunmemory_system.cpp
    sunnonlinsol/sunnonlinsol_fixedpoint.cpp
    sunnonlinsol/sunnonlinsol_newton.cpp
    test/sundials4py_test.cpp)

# Create the Python sundials library
nanobind_add_module(sundials4py STABLE_ABI NB_STATIC ${sundials_SOURCES})

# Include private header locations
target_include_directories(
  sundials4py
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include ${SUNDIALS_SOURCE_DIR}/src
          ${SUNDIALS_SOURCE_DIR}/src/arkode ${SUNDIALS_SOURCE_DIR}/src/cvodes
          ${SUNDIALS_SOURCE_DIR}/src/idas ${SUNDIALS_SOURCE_DIR}/src/kinsol)

# Link against sundials libraries
target_link_libraries(
  sundials4py
  PRIVATE sundials_arkode
          sundials_cvodes
          sundials_idas
          sundials_kinsol
          sundials_nvecserial
          sundials_sunlinsolspgmr
          sundials_sunlinsoldense
          sundials_sunlinsolband
          sundials_sunlinsolspbcgs
          sundials_sunlinsolspfgmr
          sundials_sunlinsolsptfqmr
          sundials_sunlinsolpcg
          sundials_sunmatrixband
          sundials_sunmatrixdense
          sundials_sunmatrixsparse
          sundials_sundomeigestpower
          sundials_core)

# Install directive for scikit-build-core
install(TARGETS sundials4py LIBRARY DESTINATION .)
