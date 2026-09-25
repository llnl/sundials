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
# CVODE Ginkgo heat equation example
set(cpu_gpu_examples "cv_heat_2d_ginkgo.cpp\;\;develop")
sundials_add_examples_ginkgo(cpu_gpu_examples TARGETS sundials_cvode
                              BACKENDS REF OMP CUDA HIP SYCL)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out cv_heat_2d_ginkgo*.out)
  install(FILES cv_heat_2d_ginkgo.cpp cv_heat_2d_ginkgo.hpp ${example_out}
          "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_heat_2d/cpp-ginkgo)
endif()
