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
# CVODE oneMKL C++ example
set(CVODE_examples "cvRoberts_blockdiag_onemkl.cpp\;\;develop")
sundials_add_executable(cvRoberts_blockdiag_onemkl cvRoberts_blockdiag_onemkl.cpp)
set_target_properties(cvRoberts_blockdiag_onemkl PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_blockdiag_onemkl
  PRIVATE sundials_cvode sundials_nvecsycl sundials_sunlinsolonemkldense
          MKL::MKL_DPCPP ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_blockdiag_onemkl cvRoberts_blockdiag_onemkl
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_blockdiag_onemkl.out EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_blockdiag_onemkl.cpp cvRoberts_blockdiag_onemkl.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/cpp-onemkl)
endif()
