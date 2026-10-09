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
foreach(example cv_roberts_dns cv_roberts_dns_constraints cv_roberts_dns_negsol
                cv_roberts_dns_uw)
  set(example_options)
  if(example STREQUAL "cv_roberts_dns_negsol")
    list(APPEND example_options EXAMPLE_TYPE exclude-single)
  endif()
  if(example STREQUAL "cv_roberts_dns")
    list(APPEND example_options INSTALL_FILES cv_roberts_dns_stats.csv)
  endif()
  sundials_add_example(${example} ${example}.c ${example_options})
  target_link_libraries(${example} PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()
