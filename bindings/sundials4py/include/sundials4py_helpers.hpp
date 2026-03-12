/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_HELPERS_HPP
#define _SUNDIALS4PY_HELPERS_HPP

#include "sundials4py.hpp"

namespace nb = nanobind;

namespace sundials4py {

/// This function will call a user-supplied Python function through C++ side wrappers
/// \tparam FnType is the function signature, e.g., std::remove_pointer_t<CVRhsFn>
/// \tparam FnTableType is the struct function table that holds the user-supplied Python functions as std::function
/// \tparam UserDataArg is the index of the void* user_data argument of the C function.
///         We are counting from the last arg to the first arg, so if user_data is the last arg then this should be 1.
/// \tparam Args is the template parameter pack that contains all of the types of the function arguments to the C function
///
/// \param fn_member is the name of the function in the FnTableType to call
/// \param args is the arguments to the C function, which will be forwarded to the user-supplied Python function,
///        except user_data, which is intercepted and passed as a nullptr.
template<typename FnType, typename FnTableType, std::size_t UserDataArg, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::*fn_member, Args... args)
{
  constexpr size_t N            = sizeof...(Args);
  constexpr int user_data_index = N - UserDataArg;
  auto args_tuple               = std::tuple<Args...>(args...);

  // Extract user_data from the specified position
  void* user_data = std::get<user_data_index>(args_tuple);

  // Cast user_data to FnTableType*
  auto fn_table = static_cast<FnTableType*>(user_data);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  // Pass nullptr as user_data since we do not want the user to mess with user_data (which holds our function table)
  std::get<user_data_index>(args_tuple) = nullptr;
  return std::apply([&](auto&&... call_args) { return fn(call_args...); },
                    args_tuple);
}

/// This function will call a user-supplied Python function through C++ side wrappers
/// \tparam FnType is the function signature, e.g., std::remove_pointer_t<CVRhsFn>
/// \tparam FnTableType is the struct function table that holds the user-supplied Python functions as std::function
/// \tparam MemType the type that user_data will be cast to
/// \tparam UserDataArg is the index of the void* user_data argument of the C function.
///         We are counting from the last arg to the first arg, so if user_data is the last arg then this should be 1.
/// \tparam Args is the template parameter pack that contains all of the types of the
///         function arguments to the C function
///
/// \param fn_member is the name of the function in the FnTableType to call
/// \param args is the arguments to the C function, which will be forwarded to the user-supplied
///        Python function, except user_data, which is intercepted and passed as a nullptr.
template<typename FnType, typename FnTableType, typename MemType,
         std::size_t UserDataArg, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::*fn_member, Args... args)
{
  constexpr size_t N = sizeof...(Args);
  static_assert(UserDataArg >= 1 && UserDataArg <= N);
  constexpr size_t user_data_index = N - UserDataArg;
  auto args_tuple                  = std::tuple<Args...>(args...);

  // Extract user_data from the specified position
  void* user_data = std::get<user_data_index>(args_tuple);

  // Cast user_data to FnTableType*
  auto mem      = static_cast<MemType>(user_data);
  auto fn_table = static_cast<FnTableType*>(mem->python);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  // Pass nullptr as user_data since we do not want the user to mess with
  // user_data (which holds our function table)
  static_assert(
    std::is_same_v<std::tuple_element_t<user_data_index, decltype(args_tuple)>,
                   void*>);
  std::get<user_data_index>(args_tuple) = nullptr;
  return std::apply([&](auto&&... call_args) { return fn(call_args...); },
                    args_tuple);
}

/// This function will call a user-supplied Python function through C++ side wrappers
/// \tparam FnType is the function signature, e.g., std::remove_pointer_t<CVRhsFn>
/// \tparam FnTableType is the struct function table that holds the user-supplied Python functions as std::function
/// \tparam Args is the template parameter pack that contains all of the types of the function arguments to the C function
///
/// \param fn_member is the name of the function in the FnTableType to call
/// \param args is the arguments to the C function, which will be forwarded to the user-supplied Python function.
template<typename FnType, typename FnTableType, typename T, typename... Args>
int user_supplied_fn_caller(nb::object FnTableType::*fn_member, Args... args)
{
  auto args_tuple = std::tuple<Args...>(args...);

  // Cast object->python to FnTableType*
  auto object   = static_cast<T>(std::get<0>(args_tuple));
  auto fn_table = static_cast<FnTableType*>(object->python);
  auto fn       = nb::cast<std::function<FnType>>(fn_table->*fn_member);

  return std::apply([&](auto&&... call_args) { return fn(call_args...); },
                    args_tuple);
}

///
/// \brief Helper struct to manage reference lifetimes for function return values in Python bindings.
///
/// Enables the nb::keep_alive<Nurse, Patient> paradigm when the function returns a sequence where
/// elements of the sequence (e.g., a tuple) are Nurses.
///
/// \tparam IP Index of the input argument who is kept alive by the nurses.
///         The range of IP is [0, num_args] where the 0 is the return value and 1...num_args
///         are the function input arguments in natural order.
/// \tparam IN Indexes of the return values in the returned sequence who keep the patient alive.
///         The range of IN is [0,num_tuple_entries) where the 0 is the 0th entry in the tuple.
///
/// See https://nanobind.readthedocs.io/en/latest/api_core.html#_CPPv4I0EN8nanobind11call_policyE.
///
template<size_t IP, size_t... IN>
struct returns_references_to
{
  static void precall(PyObject**, size_t, nb::detail::cleanup_list*) {}

  template<size_t N>
  static void postcall(PyObject** args, std::integral_constant<size_t, N>,
                       nb::handle ret)
  {
    static_assert(IP > 0 && IP <= N,
                  "IP in returns_references_to<IP> must be in the "
                  "range [1, number of C++ function arguments]");

    if (!nb::isinstance<nb::sequence>(ret))
    {
      throw sundials4py::error_returned("return value should be a sequence");
    }

    // Directly apply keep_alive for each IN using a fold expression
    (nb::detail::keep_alive(ret[IN].ptr(), args[IP - 1]), ...);
  }
};

} // namespace sundials4py

#endif
