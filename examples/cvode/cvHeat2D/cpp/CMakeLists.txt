# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
# CVODE serial C++ heat equation example
set(CVODE_examples "cv_heat_2d.cpp\;\;\;\;exclude-single")

sundials_add_executable(cv_heat_2d cv_heat_2d.cpp)
set_target_properties(cv_heat_2d PROPERTIES FOLDER "Examples")
if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
  target_include_directories(cv_heat_2d PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)
else()
  target_include_directories(cv_heat_2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
target_link_libraries(cv_heat_2d PRIVATE sundials_cvode sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_heat_2d PRIVATE sundials_cvode_fused_stubs)
endif()
sundials_add_test(cv_heat_2d cv_heat_2d ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_heat_2d.out EXAMPLE_TYPE exclude-single)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_heat_2d.cpp cv_heat_2d.hpp cv_heat_2d.out plot_heat_2d.py
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvHeat2D/cpp)
endif()
