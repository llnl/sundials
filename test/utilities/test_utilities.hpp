/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * Utility functions for C++ unit tests.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

// Check function return flag
inline int check_flag(const int flag, const std::string& funcname)
{
  if (!flag) return 0;
  if (flag < 0) std::cerr << "ERROR: ";
  std::cerr << funcname << " returned " << flag << std::endl;
  return 1;
}

// Check if a function returned a NULL pointer
inline int check_ptr(const void* ptr, const std::string& funcname)
{
  if (ptr) return 0;
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     sunrealtype& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    dest = stof(*(it + 1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    dest = stod(*(it + 1));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    dest = stold(*(it + 1));
#endif
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     long long& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     long int& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stol(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     int& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     bool& dest, bool store = true)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     std::string& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = *(it + 1);
    args.erase(it, it + 2);
  }
}
