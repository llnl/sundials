# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# matrix level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native matrices
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(sparse)

# optional TPL matrices
if(BUILD_SUNMATRIX_CUSPARSE)
  add_subdirectory(cusparse)
endif()

if(BUILD_SUNMATRIX_GINKGO)
  install(CODE "MESSAGE(\"\nInstall SUNMATRIX_GINKGO\n\")")
  message(STATUS "Added BUILD_SUNMATRIX_GINKGO module")
  add_library(sundials_sunmatrixginkgo INTERFACE)
  target_link_libraries(sundials_sunmatrixginkgo INTERFACE sundials_core)
  target_include_directories(
    sundials_sunmatrixginkgo
    INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
              $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
              $<BUILD_INTERFACE:${SUNDIALS_SOURCE_DIR}/src/sundials>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  install(
    FILES ${PROJECT_SOURCE_DIR}/include/sunmatrix/sunmatrix_ginkgo.hpp
          ${PROJECT_SOURCE_DIR}/include/sunmatrix/sunmatrix_ginkgobatch.hpp
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sunmatrix")
  set(_SUNDIALS_INSTALLED_COMPONENTS
      "sunmatrixginkgo;${_SUNDIALS_INSTALLED_COMPONENTS}"
      CACHE INTERNAL "" FORCE)
endif()

if(BUILD_SUNMATRIX_KOKKOSDENSE)
  install(CODE "MESSAGE(\"\nInstall SUNMATRIX_KOKKOSDENSE\n\")")
  message(STATUS "Added BUILD_SUNMATRIX_KOKKOSDENSE module")
  add_library(sundials_sunmatrixkokkosdense INTERFACE)
  target_link_libraries(
    sundials_sunmatrixkokkosdense INTERFACE sundials_core Kokkos::kokkos
                                            Kokkos::kokkoskernels)
  target_include_directories(
    sundials_sunmatrixkokkosdense
    INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
              $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
              $<BUILD_INTERFACE:${SUNDIALS_SOURCE_DIR}/src/sundials>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  install(
    FILES ${PROJECT_SOURCE_DIR}/include/sunmatrix/sunmatrix_kokkosdense.hpp
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sunmatrix")
  set(_SUNDIALS_INSTALLED_COMPONENTS
      "sunmatrixkokkosdense;${_SUNDIALS_INSTALLED_COMPONENTS}"
      CACHE INTERNAL "" FORCE)
endif()

if(BUILD_SUNMATRIX_MAGMADENSE)
  add_subdirectory(magmadense)
endif()

if(BUILD_SUNMATRIX_ONEMKLDENSE)
  add_subdirectory(onemkldense)
endif()

if(BUILD_SUNMATRIX_SLUNRLOC)
  add_subdirectory(slunrloc)
endif()
