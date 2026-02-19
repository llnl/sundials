# ---------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------

# ~~~
# sundials_option(<variable> <type> <docstring> <default value>
#                 [DEPENDS_ON dependencies]
#                 [DEPNDS_ON_THROW_ERROR])
# ~~~
#
# Within CMake creates a cache variable <variable> and sets it to the value
# <default value> if <variable> is not yet defined and, if provided, all of its
# dependencies evaluate to true. Otherwise, <variable> is not created. <type>
# may be any of the types valid for CMake's set command (FILEPATH, PATH, STRING,
# BOOL, INTERNAL). <docstring> is a description of the <variable>.
#
# The DEPENDS_ON option can be used to provide variables which must evaluate to
# true for <variable> to be created. If the dependencies do not all evaluate to
# true and <variable> exists, then a warning is printed and <variable> is unset.
#
# The DEPENDS_ON_THROW_ERROR option will change the warning to be an error.
#
# The OPTIONS option can be used to provide a list of valid <variable> values.
#
# The ADVANCED option can be used to make <variable> an advanced CMake option.

macro(sundials_option NAME TYPE DOCSTR DEFAULT_VALUE)

  # macro options and keyword inputs followed by multiple values
  set(options DEPENDS_ON_THROW_ERROR ADVANCED)
  set(multiValueArgs OPTIONS DEPENDS_ON)

  # parse inputs and create variables sundials_option_<keyword>
  cmake_parse_arguments(sundials_option "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})

  # check if dependencies for this option have been met
  set(all_depends_on_dependencies_met TRUE)
  if(sundials_option_DEPENDS_ON)
    foreach(_dependency ${sundials_option_DEPENDS_ON})
      if(NOT ${_dependency})
        set(all_depends_on_dependencies_met FALSE)
        list(APPEND depends_on_dependencies_not_met "${_dependency},")
      endif()
    endforeach()
  endif()

  if(all_depends_on_dependencies_met)

    if(NOT DEFINED ${NAME})
      set(${NAME}
          "${DEFAULT_VALUE}"
          CACHE ${TYPE} ${DOCSTR})
    else()
      set(${NAME}
          "${${NAME}}"
          CACHE ${TYPE} ${DOCSTR})
    endif()

    # make the option advanced if necessary
    if(sundials_option_ADVANCED)
      mark_as_advanced(FORCE ${NAME})
    endif()

  else()

    # if necessary, remove the CACHE variable i.e., all the variable
    # dependencies were previously met but are no longer satisfied
    if(DEFINED ${NAME})
      string(
        CONCAT
          _warn_msg_string
          "The variable ${NAME} was set to ${${NAME}} but not all of its "
          "dependencies (${depends_on_dependencies_not_met}) evaluate to TRUE. "
          "Unsetting ${NAME}.")
      unset(${NAME} CACHE)
      if(sundials_option_DEPENDS_ON_THROW_ERROR)
        message(FATAL_ERROR "${_warn_msg_string}")
      else()
        message(WARNING "${_warn_msg_string}")
      endif()
    endif()

  endif()

  # check for valid option choices
  if((DEFINED ${NAME}) AND sundials_option_OPTIONS)
    foreach(_option ${${NAME}})
      if(NOT (${_option} IN_LIST sundials_option_OPTIONS))
        list(JOIN sundials_option_OPTIONS ", " _options_msg)
        message(FATAL_ERROR "Value of ${NAME} must be one of ${_options_msg}")
      endif()
    endforeach()
    get_property(
      is_in_cache
      CACHE ${NAME}
      PROPERTY TYPE)
    if(is_in_cache)
      set_property(CACHE ${NAME} PROPERTY STRINGS ${sundials_option_OPTIONS})
    endif()
    unset(is_in_cache)
  endif()

  unset(all_depends_on_dependencies_met)
  unset(depends_on_dependencies_not_met)

endmacro()
