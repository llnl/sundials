# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

# ~~~
# sundials_install_examples_ginkgo(<MODULE>
#                                  [CPU_EXAMPLES_VAR var]
#                                  [GPU_EXAMPLES_VAR var]
#                                  [CPU_GPU_EXAMPLES_VAR var]
#                                  [DESTINATION path]
#                                  [SUNDIALS_COMPONENTS components]
#                                  [SUNDIALS_TARGETS targets]
#                                  [DEPENDENCIES files]
#                                  [EXTRA_FILES files])
# ~~~
#
# Adds an install target for each example tuple in CPU_EXAMPLES_VAR,
# GPU_EXAMPLES_VAR, and CPU_GPU_EXAMPLES_VAR that go with MODULE (e.g. cvode,
# sunlinsol).
#
# The DESTINATION option is the path *within* SUNDIALS_EXAMPLES_INSTALL_PATH
# that the files should be installed under.
#
# The SUNDIALS_COMPONENTS option is a list of CMake targets in the SUNDIALS::
# namespace provided to find_package. Note this may be the same as or a superset
# of SUNDIALS_TARGETS depending on the CMakeLists.txt template.
#
# The SUNDIALS_TARGETS option is a list of CMake targets in the SUNDIALS::
# namespace provided to target_link_libraries. Note this may be the same as or a
# subset of SUNDIALS_COMPONENTS depending on the CMakeLists.txt template.
#
# The DEPENDENCIES option is a list of additional source files that the examples
# are dependent on.
#
# The EXTRA_FILES option is a list of files to install that are not example
# source code.

macro(sundials_install_examples_ginkgo MODULE)

  set(options)
  set(oneValueArgs DESTINATION)
  set(multiValueArgs
      CPU_EXAMPLES_VAR
      GPU_EXAMPLES_VAR
      CPU_GPU_EXAMPLES_VAR
      SUNDIALS_COMPONENTS
      SUNDIALS_TARGETS
      EXTRA_FILES
      DEPENDENCIES)

  # Parse keyword arguments/options
  cmake_parse_arguments(arg "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # Install the example source, header, and output file
  foreach(example_type CPU GPU CPU_GPU)
    foreach(example_tuple ${${arg_${example_type}_EXAMPLES_VAR}})

      # filename must be the first item in the list of example tuples
      list(GET example_tuple 0 example)

      # extract the file name without extension
      get_filename_component(example_noext ${example} NAME_WE)

      # get example header and output files
      file(GLOB example_header ${example_noext}.h*)
      file(GLOB example_out ${example_noext}*.out)

      # install files
      install(FILES ${example} ${example_header} ${example_out}
              DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/${arg_DESTINATION})

    endforeach()
  endforeach()

  # Install the extra files and dependencies
  if(arg_EXTRA_FILES OR arg_DEPENDENCIES)
    install(FILES ${arg_EXTRA_FILES} ${arg_DEPENDENCIES}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/${arg_DESTINATION})
  endif()

  # Prepare substitution variables for CMakeLists and/or Makefile templates
  if(arg_DEPENDENCIES)
    list2string(arg_DEPENDENCIES EXAMPLES_DEPENDENCIES)
  endif()

  if(arg_CPU_EXAMPLES_VAR)
    examples2string(${arg_CPU_EXAMPLES_VAR} CPU_EXAMPLES)
  endif()
  if(arg_GPU_EXAMPLES_VAR)
    examples2string(${arg_GPU_EXAMPLES_VAR} GPU_EXAMPLES)
  endif()
  if(arg_CPU_GPU_EXAMPLES_VAR)
    examples2string(${arg_CPU_GPU_EXAMPLES_VAR} CPU_GPU_EXAMPLES)
  endif()

  # components for find_package
  if(arg_SUNDIALS_COMPONENTS)
    list2string(arg_SUNDIALS_COMPONENTS EXAMPLES_CMAKE_COMPONENTS)
  endif()

  # targets for target_link_libraries
  foreach(target ${arg_SUNDIALS_TARGETS})
    list(APPEND target_list SUNDIALS::${target})
  endforeach()
  if(target_list)
    list2string(target_list EXAMPLES_CMAKE_TARGETS)
  endif()

  # Generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_CXX_ginkgo_ex.in
    ${PROJECT_BINARY_DIR}/examples/${arg_DESTINATION}/CMakeLists.txt @ONLY)

  # Install CMakelists.txt
  install(FILES ${PROJECT_BINARY_DIR}/examples/${arg_DESTINATION}/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/${arg_DESTINATION})

  # Add test_install target
  sundials_add_test_install(${MODULE} ginkgo)

endmacro()
