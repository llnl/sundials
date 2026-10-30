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
# Generate a standalone CMakeLists.txt for an installed example directory.
# ---------------------------------------------------------------

function(sundials_install_example_cmake _target)
  set(source_dir "${CMAKE_CURRENT_SOURCE_DIR}")
  file(RELATIVE_PATH _example_subdir "${PROJECT_SOURCE_DIR}/examples"
       "${source_dir}")
  string(REGEX MATCH "^[^/]+" package "${_example_subdir}")
  set(destination "${SUNDIALS_EXAMPLES_INSTALL_PATH}/${_example_subdir}")

  get_property(_languages DIRECTORY PROPERTY SUNDIALS_EXAMPLE_LANGUAGES)
  get_property(_target_commands DIRECTORY PROPERTY SUNDIALS_EXAMPLE_COMMANDS)
  set(_compiler_settings)
    get_target_property(_sources ${_target} SOURCES)
    get_target_property(_target_source_dir ${_target} SOURCE_DIR)
    set(_installed_sources)
    set(_source_commands)
    foreach(_source IN LISTS _sources)
      if(_source MATCHES "^\\$<")
        list(APPEND _installed_sources "${_source}")
        continue()
      endif()
      get_filename_component(_installed_source "${_source}" NAME)
      list(APPEND _installed_sources "${_installed_source}")

      # Preserve the target's source spelling for scoped property lookup, and
      # fall back to its absolute path when the relative path is interpreted
      # from this function's directory.
      if(NOT IS_ABSOLUTE "${_source}")
        set(_absolute_source "${_target_source_dir}/${_source}")
      endif()

      get_source_file_property(_language "${_source}" TARGET_DIRECTORY
                               ${_target} LANGUAGE)
      if(_language STREQUAL "NOTFOUND" AND NOT IS_ABSOLUTE "${_source}")
        get_source_file_property(_language "${_absolute_source}"
                                 TARGET_DIRECTORY ${_target} LANGUAGE)
      endif()
      if(_language)
        list(APPEND _languages "${_language}")
        get_filename_component(_extension "${_source}" LAST_EXT)
        if(_extension STREQUAL ".cpp" AND NOT _language STREQUAL "CXX")
          string(
            APPEND
            _source_commands
            "set_source_files_properties(${_installed_source} PROPERTIES LANGUAGE ${_language})\n"
          )
        endif()
      elseif(_source MATCHES "\\.cu$")
        list(APPEND _languages CUDA)
      elseif(_source MATCHES "\\.(cc|cpp|cxx)$")
        list(APPEND _languages CXX)
      elseif(_source MATCHES "\\.(f|f90|F|F90)$")
        list(APPEND _languages Fortran)
      elseif(_source MATCHES "\\.c$")
        list(APPEND _languages C)
      endif()

    endforeach()

    string(APPEND _target_commands "${_source_commands}")
    string(JOIN " " _source_list ${_installed_sources})
    string(APPEND _target_commands
           "add_executable(${_target} ${_source_list})\n")

    get_target_property(_include_dirs ${_target} INCLUDE_DIRECTORIES)
    if(_include_dirs)
      foreach(_include_dir IN LISTS _include_dirs)
        if(_include_dir MATCHES "/examples/utilities$")
          string(
            APPEND _target_commands
            "target_include_directories(${_target} PRIVATE \${CMAKE_CURRENT_SOURCE_DIR})\n"
          )
        endif()
      endforeach()
    endif()

    get_target_property(_definitions ${_target} COMPILE_DEFINITIONS)
    if(_definitions)
      set(_definition_list)
      foreach(_definition IN LISTS _definitions)
        list(APPEND _definition_list "[==[${_definition}]==]")
      endforeach()
      string(JOIN " " _definition_list ${_definition_list})
      string(
        APPEND _target_commands
        "target_compile_definitions(${_target} PRIVATE ${_definition_list})\n")
    endif()

    get_target_property(_link_libraries ${_target} LINK_LIBRARIES)
    if(_link_libraries)
      set(_installed_link_libraries)
      foreach(_library IN LISTS _link_libraries)
        if(_library MATCHES "^MPI::MPI_(C|CXX|Fortran)$")
          list(APPEND _languages "${CMAKE_MATCH_1}")
        endif()
        if(_library MATCHES "^sundials_(.+)$")
          list(APPEND _installed_link_libraries "SUNDIALS::${CMAKE_MATCH_1}")
        else()
          list(APPEND _installed_link_libraries "${_library}")
        endif()
      endforeach()
      set(_library_list)
      foreach(_library IN LISTS _installed_link_libraries)
        list(APPEND _library_list "[==[${_library}]==]")
      endforeach()
      string(JOIN " " _library_list ${_library_list})
      string(APPEND _target_commands
             "target_link_libraries(${_target} PRIVATE ${_library_list})\n")
    endif()

    foreach(_property IN ITEMS CUDA_ARCHITECTURES CXX_STANDARD LINKER_LANGUAGE)
      get_target_property(_value ${_target} ${_property})
      if(_value)
        if(_property STREQUAL "CXX_STANDARD")
          # Accelerator sources may consume C++ compile features from
          # dependencies. Ignore the globally initialized C++ standard for
          # targets that do not use C++ or an accelerator language.
          if(CUDA IN_LIST _languages OR HIP IN_LIST _languages)
            list(APPEND _languages CXX)
          elseif(NOT CXX IN_LIST _languages)
            continue()
          endif()
        endif()
        string(
          APPEND _target_commands
          "set_property(TARGET ${_target} PROPERTY ${_property} ${_value})\n")
      endif()
    endforeach()
    if(Fortran IN_LIST _languages)
      string(
        APPEND
        _target_commands
        "set_property(TARGET ${_target} PROPERTY Fortran_MODULE_DIRECTORY \${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${_target}.dir)\n"
      )
    endif()
    string(APPEND _target_commands
           "add_test(NAME ${_target} COMMAND ${_target})\n\n")
  list(REMOVE_DUPLICATES _languages)
  set_property(DIRECTORY PROPERTY SUNDIALS_EXAMPLE_LANGUAGES "${_languages}")
  set_property(DIRECTORY PROPERTY SUNDIALS_EXAMPLE_COMMANDS "${_target_commands}")
  foreach(_language IN LISTS _languages)
    set(_variable CMAKE_${_language}_COMPILER)
    if(DEFINED ${_variable} AND NOT "${${_variable}}" STREQUAL "")
      string(
        APPEND
        _compiler_settings
        "set(${_variable}\n  [==[${${_variable}}]==]\n  CACHE FILEPATH \"${_language} compiler\")\n\n"
      )
    endif()

    set(_variable CMAKE_${_language}_FLAGS)
    if(DEFINED ${_variable})
      string(
        APPEND
        _compiler_settings
        "set(${_variable}\n  [==[${${_variable}}]==]\n  CACHE STRING \"${_language} compiler flags\")\n\n"
      )
    endif()

    set(_variable CMAKE_${_language}_STANDARD)
    if(DEFINED ${_variable} AND NOT "${${_variable}}" STREQUAL "")
      string(
        APPEND
        _compiler_settings
        "set(${_variable}\n  [==[${${_variable}}]==]\n  CACHE STRING \"${_language} standard\")\n\n"
      )
    endif()
  endforeach()

  if(CUDA IN_LIST _languages
     AND DEFINED CMAKE_CUDA_HOST_COMPILER
     AND NOT "${CMAKE_CUDA_HOST_COMPILER}" STREQUAL "")
    string(
      APPEND
      _compiler_settings
      "set(CMAKE_CUDA_HOST_COMPILER\n  [==[${CMAKE_CUDA_HOST_COMPILER}]==]\n  CACHE FILEPATH \"CUDA host compiler\")\n\n"
    )
  endif()

  string(JOIN " " _language_list ${_languages})
  set(PACKAGE "${package}")
  set(COMPILER_SETTINGS "${_compiler_settings}")
  set(LANGUAGES "${_language_list}")
  set(SUNDIALS_CMAKE_DIR
      "${CMAKE_INSTALL_PREFIX}/${SUNDIALS_INSTALL_CMAKEDIR}")
  if(Fortran IN_LIST _languages)
    string(PREPEND _target_commands "set(CMAKE_Fortran_PREPROCESS ON)\n\n")
  endif()
  set(TARGET_COMMANDS "${_target_commands}")

  set(_output "${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt")
  configure_file(
    "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/SundialsInstallExampleCMakeLists.txt.in"
    "${_output}" @ONLY)
  install(FILES "${_output}" DESTINATION "${destination}")
endfunction()
