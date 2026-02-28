/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * ---------------------------------------------------------------------------*/

#include <chrono>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <string>
#include <thread>

#include "sundials/sundials_math.h"
#include "sundials/sundials_profiler.h"
#include "sundials/sundials_types.h"

static int sleep(SUNProfiler prof, int sec, double* chrono)
{
  auto begin = std::chrono::steady_clock::now();

  // We dont check the flag returned to avoid introducing extra
  // overhead which makes it harder to compare the SUNProfiler
  // time with the chrono time.
  SUNProfiler_Begin(prof, "sleep");
  std::this_thread::sleep_for(std::chrono::seconds(sec));
  SUNProfiler_End(prof, "sleep");

  auto end = std::chrono::steady_clock::now();

  auto elapsed =
    std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
  *chrono = std::chrono::duration<double>(elapsed).count() * 1e-9;

  return 0;
}

static int print_timings(SUNProfiler prof)
{
  // Output timing in default (table) format
  int flag = SUNProfiler_Print(prof, stdout);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_Print returned " << flag << "\n";
    return 1;
  }

  return 0;
}

int main()
{
  // -----
  // Setup
  // -----

  double chrono;

  std::cout << "Testing SUNProfiler\n";

  SUNProfiler prof = nullptr;
  int flag = SUNProfiler_Create(SUN_COMM_NULL, "SUNProfiler Test", &prof);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_Create returned " << flag << "\n";
    return 1;
  }

  // ------
  // Test 1
  // ------

  std::cout << "\nTest 1: sleep 1s, check timings \n";

  flag = sleep(prof, 1, &chrono);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "sleep returned " << flag << "\n";
    return 1;
  }

  flag = print_timings(prof);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "print_timings returned " << flag << "\n";
    return 1;
  }

  double time       = 0;
  double resolution = 0;
  flag              = SUNProfiler_GetElapsedTime(prof, "sleep", &time);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_GetElapsedTime returned " << flag << "\n";
    return 1;
  }

  flag = SUNProfiler_GetTimerResolution(prof, &resolution);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_GetTimerResolution returned " << flag << "\n";
    return 1;
  }

  if (SUNRCompareTol(time, chrono, 1e-2))
  {
    std::cerr << ">>> FAILURE: "
              << "time recorded was " << time << "s, but expected " << chrono
              << "s +/- " << 1e-2 << "\n";
    return 1;
  }

  // ------
  // Test 2
  // ------

  std::cout << "\nTest 2: reset, sleep 2s, print timings\n";

  // Reset timing
  flag = SUNProfiler_Reset(prof);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_Reset returned " << flag << "\n";
    return 1;
  }

  flag = sleep(prof, 2, &chrono);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "sleep returned " << flag << "\n";
    return 1;
  }

  flag = print_timings(prof);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "print_timings returned " << flag << "\n";
    return 1;
  }

  flag = SUNProfiler_GetElapsedTime(prof, "sleep", &time);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_GetElapsedTime returned " << flag << "\n";
    return 1;
  }

  flag = SUNProfiler_GetTimerResolution(prof, &resolution);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_GetTimerResolution returned " << flag << "\n";
    return 1;
  }

  if (SUNRCompareTol(time, chrono, 1e-2))
  {
    std::cerr << ">>> FAILURE: "
              << "time recorded was " << time << "s, but expected " << chrono
              << "s +/- " << 1e-2 << "\n";
    return 1;
  }

  // ------
  // Test 3
  // ------

  std::cout << "\nTest 3: multiple outputs to a file for plot test\n";

  std::FILE* fout = std::fopen("profiling_test_output.txt", "w");
  if (fout == nullptr)
  {
    std::cerr << ">>> FAILURE: "
              << "fopen returned a null pointer\n";
    return 1;
  }

  for (int i = 1; i < 4; i++)
  {
    // Reset timing
    flag = SUNProfiler_Reset(prof);
    if (flag)
    {
      std::cerr << ">>> FAILURE: "
                << "SUNProfiler_Reset returned " << flag << "\n";
      return 1;
    }

    flag = sleep(prof, i, &chrono);
    if (flag)
    {
      std::cerr << ">>> FAILURE: "
                << "sleep returned " << flag << "\n";
      return 1;
    }

    flag = SUNProfiler_Print(prof, fout);
    if (flag)
    {
      std::cerr << ">>> FAILURE: "
                << "SUNProfiler_Print returned " << flag << "\n";
      return 1;
    }
  }

  std::fclose(fout);

  // --------
  // Clean up
  // --------

  flag = SUNProfiler_Free(&prof);
  if (flag)
  {
    std::cerr << ">>> FAILURE: "
              << "SUNProfiler_Free returned " << flag << "\n";
    return 1;
  }

  std::cout << "\nSUCCESS - test complete\n";

  return 0;
}
