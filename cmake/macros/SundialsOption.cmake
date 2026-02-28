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

#[=======================================================================[.rst:
SundialsOption
--------------

This module provides a command for setting SUNDIALS cache variables
(configuration options).

Load this module in with:

.. code-block:: cmake

   include(SundialsOption)

Commands
^^^^^^^^

This module provides the following command:

.. cmake:command:: sundials_option

   Set a SUNDIALS cache variable (configuration option) to a given value.

   .. code-block:: cmake

      sundials_option(<variable> <type> <help string> <default value>
                      [ADVANCED]
                      [OPTIONS options]
                      [DEPENDS_ON dependencies]
                      [DEPENDS_ON_THROW_ERROR]
                      [DEPRECATED_NAMES names]
                      [NEGATE_DEPRECATED])

   Wraps the CMake :cmake:command:`set() <cmake:command:set>` command to set the
   given cache variable.

   The arguments are:

   ``<variable>``
     The name of a variable that stores the option value.

   ``<type>``
     The type of the cache entry.

   ``<help string>``
     Text providing a quick summary of the option for CMake GUIs.

   ``<default value>``
     The default value for the cache variable.

     * If the cache entry does not exist prior to the call, then it will be set
       to the default value.

   The options are:

   ``ADVANCED``
     Mark the cache variable as advanced.

   ``OPTIONS <options>...``
     A list of valid values for the cache variable.

     * If a user set the variable to a value not in the list of options, the
       configuration is halted and an error message is printed.

   ``DEPENDS_ON <dependencies>...``
     A list of variables which must evaluate to true for the cache variable to
     be set.

     * If any of the dependencies evaluate to false, then the cache variable will
       be unset and a warning message is printed.

   ``DEPENDS_ON_THROW_ERROR``
     Throw an error if the option dependencies are not met.

   ``DEPRECATED_NAMES <variables>...``
     A list of deprecated cache variable names for the cache variable.

     * If the cache variable is already defined and any of the deprecated
       cache variables are also defined, the deprecated variables are
       ignored and a warning message is printed.

     * If the cache variable is not defined and a deprecated cache variable is
       defined, the value of the deprecated variable is copied to the cache
       variable and a warning message is printed.

     * If multiple deprecated cache variable names are provided and more than
       one of them is defined, an error is printed if the deprecated variables
       do not all have the same value.

     * If ``SUNDIALS_ENABLE_UNSET_DEPRECATED`` is true, any deprecated variables
       will be unset after checking the value for the cases above.

   ``NEGATE_DEPRECATED``
     Negate the value of deprecated boolean cache variables when setting a cache
     variable that is not already defined.
#]=======================================================================]

function(sundials_option NAME TYPE DOCSTR DEFAULT_VALUE)

  # macro options and keyword inputs followed by multiple values
  set(options DEPENDS_ON_THROW_ERROR ADVANCED NEGATE_DEPRECATED)
  set(multiValueArgs OPTIONS DEPENDS_ON DEPRECATED_NAMES)

  # parse inputs and create variables arg_<keyword>
  cmake_parse_arguments(arg "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # check for deprecated options
  if(arg_DEPRECATED_NAMES)
    unset(_save_name)
    foreach(_deprecated_name ${arg_DEPRECATED_NAMES})
      if(DEFINED ${_deprecated_name})
        message(
          WARNING
            "The option ${_deprecated_name} is deprecated. Use ${NAME} instead."
        )
        if(NOT DEFINED _save_name)
          # save name and value separately in case unset below
          set(_save_name ${_deprecated_name})
          set(_save_value ${${_deprecated_name}})
        else()
          if(TYPE STREQUAL BOOL)
            # Check if boolean values match
            if(NOT (${_save_value} AND ${${_deprecated_name}})
               AND (${_save_value} OR ${${_deprecated_name}}))
              message(
                FATAL_ERROR
                  "Inconsistent deprecated options: ${_save_name} = ${_save_value} and ${_deprecated_name} = ${${_deprecated_name}}."
              )
            endif()
          else()
            # Check if filepath/string/path match
            if(NOT (${_save_value} STREQUAL ${${_deprecated_name}}))
              message(
                FATAL_ERROR
                  "Inconsistent deprecated options: ${_save_name} = ${_save_value} and ${_deprecated_name} = ${${_deprecated_name}}."
              )
            endif()
          endif()
          message(WARNING "Multiple deprecated options for ${NAME} provided.")
        endif()
        if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
          unset(${_deprecated_name} CACHE)
        endif()
      endif()
    endforeach()
    if(DEFINED _save_name)
      if(DEFINED ${NAME})
        message(
          WARNING
            "Both ${NAME} and ${_save_name} (deprecated) are defined. Ignoring "
            "${_save_name}.")
      else()
        if((TYPE STREQUAL BOOL) AND arg_NEGATE_DEPRECATED)
          if(${_save_value})
            set(${NAME}
                FALSE
                CACHE ${TYPE} ${DOCSTR})
          else()
            set(${NAME}
                TRUE
                CACHE ${TYPE} ${DOCSTR})
          endif()
        else()
          set(${NAME}
              "${_save_value}"
              CACHE ${TYPE} ${DOCSTR})
        endif()
      endif()
    endif()
  endif()

  # check if dependencies for this option have been met
  set(all_depends_on_dependencies_met TRUE)
  if(arg_DEPENDS_ON)
    foreach(_dependency ${arg_DEPENDS_ON})
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
          CACHE ${TYPE} ${DOCSTR} FORCE)
    endif()

    # make the option advanced if necessary
    if(arg_ADVANCED)
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
      if(arg_DEPENDS_ON_THROW_ERROR)
        message(FATAL_ERROR "${_warn_msg_string}")
      else()
        message(WARNING "${_warn_msg_string}")
      endif()
    endif()

  endif()

  # check for valid option choices
  if((DEFINED ${NAME}) AND arg_OPTIONS)
    foreach(_option ${${NAME}})
      if(NOT (${_option} IN_LIST arg_OPTIONS))
        list(JOIN arg_OPTIONS ", " _options_msg)
        message(FATAL_ERROR "Value of ${NAME} must be one of ${_options_msg}")
      endif()
    endforeach()
    get_property(
      is_in_cache
      CACHE ${NAME}
      PROPERTY TYPE)
    if(is_in_cache)
      set_property(CACHE ${NAME} PROPERTY STRINGS ${arg_OPTIONS})
    endif()
    unset(is_in_cache)
  endif()

  unset(all_depends_on_dependencies_met)
  unset(depends_on_dependencies_not_met)

endfunction()
