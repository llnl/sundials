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

function(sundials_install_example_cmake package source_dir destination)
  get_property(
    _targets
    DIRECTORY "${source_dir}"
    PROPERTY BUILDSYSTEM_TARGETS)

  set(_contents
      [=[# -----------------------------------------------------------------
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

cmake_minimum_required(VERSION 3.18)
@COMPILER_SETTINGS@
project(@PACKAGE@_example LANGUAGES @LANGUAGES@)

include(CTest)

@DEPENDENCY_SETTINGS@

set(SUNDIALS_DIR "@SUNDIALS_CMAKE_DIR@"
    CACHE PATH "Location of SUNDIALSConfig.cmake")
find_package(SUNDIALS REQUIRED NO_DEFAULT_PATH)

]=])

  set(_languages)
  set(_compiler_settings)
  set(_dependency_settings)
  set(_requires_hip FALSE)
  set(_requires_mpi FALSE)
  set(_target_commands)
  foreach(_target IN LISTS _targets)
    get_target_property(_type ${_target} TYPE)
    if(NOT _type STREQUAL "EXECUTABLE")
      continue()
    endif()

    get_target_property(_sources ${_target} SOURCES)
    set(_installed_sources)
    set(_source_commands)
    foreach(_source IN LISTS _sources)
      if(_source MATCHES "^\\$<")
        list(APPEND _installed_sources "${_source}")
        continue()
      endif()
      get_filename_component(_installed_source "${_source}" NAME)
      list(APPEND _installed_sources "${_installed_source}")

      get_source_file_property(_language "${_source}" TARGET_DIRECTORY
                               ${_target} LANGUAGE)
      if(_language AND NOT _language STREQUAL "NOTFOUND")
        list(APPEND _languages "${_language}")
        get_filename_component(_extension "${_source}" LAST_EXT)
        if((_extension STREQUAL ".cpp" OR _extension STREQUAL ".c")
           AND NOT _language STREQUAL "CXX"
           AND NOT _language STREQUAL "C")
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

      get_source_file_property(_cuda_format "${_source}" TARGET_DIRECTORY
                               ${_target} CUDA_SOURCE_PROPERTY_FORMAT)
      if(_cuda_format AND NOT _cuda_format STREQUAL "NOTFOUND")
        string(
          APPEND
          _source_commands
          "set_source_files_properties(${_installed_source} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT ${_cuda_format})\n"
        )
      endif()
    endforeach()

    string(APPEND _target_commands "${_source_commands}")
    string(JOIN " " _source_list ${_installed_sources})
    string(APPEND _target_commands
           "add_executable(${_target} ${_source_list})\n")

    get_target_property(_include_dirs ${_target} INCLUDE_DIRECTORIES)
    if(_include_dirs AND NOT _include_dirs STREQUAL "_include_dirs-NOTFOUND")
      set(_installed_include_dirs)
      foreach(_include_dir IN LISTS _include_dirs)
        if(_include_dir MATCHES "/examples/utilities$")
          list(APPEND _installed_include_dirs "\${CMAKE_CURRENT_SOURCE_DIR}")
        elseif(NOT _include_dir STREQUAL "${source_dir}")
          if(IS_ABSOLUTE "${_include_dir}")
            file(RELATIVE_PATH _source_relative "${PROJECT_SOURCE_DIR}"
                 "${_include_dir}")
            file(RELATIVE_PATH _binary_relative "${PROJECT_BINARY_DIR}"
                 "${_include_dir}")
            if(NOT _source_relative MATCHES "^\\.\\."
               OR NOT _binary_relative MATCHES "^\\.\\.")
              message(
                WARNING
                  "Skipping build-tree include directory '${_include_dir}' when generating the installed ${package} example at '${destination}'."
              )
              continue()
            endif()
          endif()
          list(APPEND _installed_include_dirs "${_include_dir}")
        endif()
      endforeach()
      if(_installed_include_dirs)
        set(_include_list)
        foreach(_include_dir IN LISTS _installed_include_dirs)
          string(REPLACE "\\" "\\\\" _escaped_include "${_include_dir}")
          string(REPLACE "\"" "\\\"" _escaped_include
                         "${_escaped_include}")
          list(APPEND _include_list "\"${_escaped_include}\"")
        endforeach()
        string(JOIN " " _include_list ${_include_list})
        string(
          APPEND _target_commands
          "target_include_directories(${_target} PRIVATE ${_include_list})\n")
      endif()
    endif()

    get_target_property(_definitions ${_target} COMPILE_DEFINITIONS)
    if(_definitions AND NOT _definitions STREQUAL "_definitions-NOTFOUND")
      set(_definition_list)
      foreach(_definition IN LISTS _definitions)
        string(REPLACE "\\" "\\\\" _escaped_definition "${_definition}")
        string(REPLACE "\"" "\\\"" _escaped_definition
                       "${_escaped_definition}")
        list(APPEND _definition_list "\"${_escaped_definition}\"")
      endforeach()
      string(JOIN " " _definition_list ${_definition_list})
      string(
        APPEND _target_commands
        "target_compile_definitions(${_target} PRIVATE ${_definition_list})\n")
    endif()

    get_target_property(_link_libraries ${_target} LINK_LIBRARIES)
    if(_link_libraries AND NOT _link_libraries STREQUAL
                           "_link_libraries-NOTFOUND")
      set(_installed_link_libraries)
      foreach(_library IN LISTS _link_libraries)
        if(_library MATCHES "^hip::")
          set(_requires_hip TRUE)
        elseif(_library MATCHES "^MPI::")
          set(_requires_mpi TRUE)
        endif()
        if(_library MATCHES "^sundials_(.+)$")
          list(APPEND _installed_link_libraries "SUNDIALS::${CMAKE_MATCH_1}")
        else()
          list(APPEND _installed_link_libraries "${_library}")
        endif()
      endforeach()
      set(_library_list)
      foreach(_library IN LISTS _installed_link_libraries)
        string(REPLACE "\\" "\\\\" _escaped_library "${_library}")
        string(REPLACE "\"" "\\\"" _escaped_library
                       "${_escaped_library}")
        list(APPEND _library_list "\"${_escaped_library}\"")
      endforeach()
      string(JOIN " " _library_list ${_library_list})
      string(APPEND _target_commands
             "target_link_libraries(${_target} PRIVATE ${_library_list})\n")
    endif()

    foreach(_property IN ITEMS CUDA_ARCHITECTURES CXX_STANDARD LINKER_LANGUAGE)
      get_target_property(_value ${_target} ${_property})
      if(_value AND NOT _value STREQUAL "_value-NOTFOUND")
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
  endforeach()

  list(REMOVE_DUPLICATES _languages)
  foreach(_language IN LISTS _languages)
    set(_variable CMAKE_${_language}_COMPILER)
    if(DEFINED ${_variable} AND NOT "${${_variable}}" STREQUAL "")
      set(_value "${${_variable}}")
      string(REPLACE "\\" "\\\\" _value "${_value}")
      string(REPLACE "\"" "\\\"" _value "${_value}")
      string(APPEND _compiler_settings
             "set(${_variable}\n  \"${_value}\"\n  CACHE FILEPATH \"${_language} compiler\")\n\n")
    endif()

    set(_variable CMAKE_${_language}_FLAGS)
    if(DEFINED ${_variable})
      set(_value "${${_variable}}")
      string(REPLACE "\\" "\\\\" _value "${_value}")
      string(REPLACE "\"" "\\\"" _value "${_value}")
      string(APPEND _compiler_settings
             "set(${_variable}\n  \"${_value}\"\n  CACHE STRING \"${_language} compiler flags\")\n\n")
    endif()

    set(_variable CMAKE_${_language}_STANDARD)
    if(DEFINED ${_variable} AND NOT "${${_variable}}" STREQUAL "")
      set(_value "${${_variable}}")
      string(REPLACE "\\" "\\\\" _value "${_value}")
      string(REPLACE "\"" "\\\"" _value "${_value}")
      string(APPEND _compiler_settings
             "set(${_variable}\n  \"${_value}\"\n  CACHE STRING \"${_language} standard\")\n\n")
    endif()
  endforeach()

  if(CUDA IN_LIST _languages AND DEFINED CMAKE_CUDA_HOST_COMPILER
     AND NOT "${CMAKE_CUDA_HOST_COMPILER}" STREQUAL "")
    set(_value "${CMAKE_CUDA_HOST_COMPILER}")
    string(REPLACE "\\" "\\\\" _value "${_value}")
    string(REPLACE "\"" "\\\"" _value "${_value}")
    string(APPEND _compiler_settings
           "set(CMAKE_CUDA_HOST_COMPILER\n  \"${_value}\"\n  CACHE FILEPATH \"CUDA host compiler\")\n\n")
  endif()

  if(_requires_mpi)
    string(APPEND _dependency_settings "find_package(MPI REQUIRED)\n\n")
    foreach(_mpi_language C CXX Fortran)
      set(_variable MPI_${_mpi_language}_COMPILER)
      if(DEFINED ${_variable} AND NOT "${${_variable}}" STREQUAL "")
        set(_value "${${_variable}}")
        string(REPLACE "\\" "\\\\" _value "${_value}")
        string(REPLACE "\"" "\\\"" _value "${_value}")
        string(APPEND _compiler_settings
               "set(${_variable}\n  \"${_value}\"\n  CACHE FILEPATH \"MPI ${_mpi_language} compiler\")\n\n")
      endif()
    endforeach()
  endif()

  if(_requires_hip)
    string(APPEND _dependency_settings "find_package(HIP REQUIRED)\n\n")
  endif()

  string(JOIN " " _language_list ${_languages})
  string(REPLACE "@PACKAGE@" "${package}" _contents "${_contents}")
  string(REPLACE "@COMPILER_SETTINGS@" "${_compiler_settings}" _contents
                 "${_contents}")
  string(REPLACE "@DEPENDENCY_SETTINGS@" "${_dependency_settings}" _contents
                 "${_contents}")
  string(REPLACE "@LANGUAGES@" "${_language_list}" _contents "${_contents}")
  string(
    REPLACE "@SUNDIALS_CMAKE_DIR@"
            "${CMAKE_INSTALL_PREFIX}/${SUNDIALS_INSTALL_CMAKEDIR}"
            _contents "${_contents}")
  if(Fortran IN_LIST _languages)
    string(APPEND _contents "set(CMAKE_Fortran_PREPROCESS ON)\n\n")
  endif()
  string(APPEND _contents "${_target_commands}")

  file(RELATIVE_PATH _relative_dir "${CMAKE_CURRENT_SOURCE_DIR}"
       "${source_dir}")
  set(_output "${CMAKE_CURRENT_BINARY_DIR}/${_relative_dir}/CMakeLists.txt")
  file(WRITE "${_output}" "${_contents}")
  install(FILES "${_output}" DESTINATION "${destination}")
endfunction()
