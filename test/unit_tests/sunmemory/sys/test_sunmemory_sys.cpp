/*------------------------------------------------------------------
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
 *-----------------------------------------------------------------*/

#include <iostream>
#include <sundials/sundials_core.hpp>
#include <sunmemory/sunmemory_system.h>

static int test_instance(SUNMemoryHelper helper, SUNMemoryType mem_type,
                         bool print_test_status)
{
  // Try and allocate some memory
  int N                 = 8;
  size_t bytes_to_alloc = N * sizeof(sunrealtype);
  SUNMemory some_memory = nullptr;

  if (print_test_status) { std::cout << "  SUNMemoryHelper_Alloc... \n"; }
  int retval = SUNMemoryHelper_Alloc(helper, &some_memory, bytes_to_alloc,
                                     mem_type, nullptr);
  if (retval)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_Alloc... FAILED\n";
    }
    return -1;
  }

  // Write to the memory
  sunrealtype* some_arr = static_cast<sunrealtype*>(some_memory->ptr);
  for (int i = 0; i < N; i++) { some_arr[i] = i * sunrealtype{1.0}; }
  if (print_test_status) { std::cout << "  SUNMemoryHelper_Alloc... PASSED\n"; }

  // Try and copy the memory
  if (print_test_status) { std::cout << "  SUNMemoryHelper_Copy... \n"; }
  SUNMemory other_memory = nullptr;
  SUNMemoryHelper_Alloc(helper, &other_memory, bytes_to_alloc, mem_type, nullptr);
  retval = SUNMemoryHelper_Copy(helper, other_memory, some_memory,
                                bytes_to_alloc, nullptr);
  if (retval)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_Copy... FAILED\n";
    }
    return -1;
  }
  else
  {
    sunrealtype* other_arr = static_cast<sunrealtype*>(other_memory->ptr);
    for (int i = 0; i < N; i++)
    {
      if (some_arr[i] != other_arr[i])
      {
        if (print_test_status)
        {
          std::cout << "  SUNMemoryHelper_Copy... FAILED\n";
        }
        return -1;
      }
    }
  }
  if (print_test_status) { std::cout << "  SUNMemoryHelper_Copy... PASSED\n"; }

  // Try and deallocate
  if (print_test_status) { std::cout << "  SUNMemoryHelper_Dealloc... \n"; }
  retval = SUNMemoryHelper_Dealloc(helper, some_memory, nullptr);
  if (retval)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_Dealloc... FAILED\n";
    }
    return -1;
  }
  retval = SUNMemoryHelper_Dealloc(helper, other_memory, nullptr);
  if (retval)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_Dealloc... FAILED\n";
    }
    return -1;
  }
  if (print_test_status)
  {
    std::cout << "  SUNMemoryHelper_Dealloc... PASSED\n";
  }

  // Check alloc stats
  if (print_test_status)
  {
    std::cout << "  SUNMemoryHelper_GetAllocStats... \n";
  }
  unsigned long num_allocations, num_deallocations;
  size_t bytes_allocated, bytes_high_watermark;

  retval = SUNMemoryHelper_GetAllocStats(helper, mem_type, &num_allocations,
                                         &num_deallocations, &bytes_allocated,
                                         &bytes_high_watermark);
  if (retval)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    }
    return -1;
  }
  if (print_test_status)
  {
    std::cout << "\tnum_allocations = " << num_allocations
              << " num_deallocations = " << num_deallocations
              << " bytes_allocated = " << bytes_allocated
              << " bytes_high_watermark = " << bytes_high_watermark << "\n";
  }
  if (num_allocations != 2)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    }
    if (print_test_status) { std::cout << "    num_allocations != 2\n"; }
    return -1;
  }
  if (num_deallocations != 2)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    }
    if (print_test_status) { std::cout << "    num_deallocations != 2\n"; }
    return -1;
  }
  if (bytes_allocated != 0)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    }
    if (print_test_status) { std::cout << "    bytes_allocated != 0\n"; }
    return -1;
  }
  if (bytes_high_watermark != bytes_to_alloc * 2)
  {
    if (print_test_status)
    {
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    }
    if (print_test_status) { std::cout << "    bytes_high_watermark != 0\n"; }
    return -1;
  }
  if (print_test_status)
  {
    std::cout << "  SUNMemoryHelper_GetAllocStats... PASSED\n";
  }
  return retval;
}

int main(int argc, char* argv[])
{
  sundials::Context sunctx;

  std::cout << "Testing the SUNMemoryHelper_Sys module... \n";

  std::cout << "  SUNMemoryHelper_Sys... \n";
  SUNMemoryHelper helper = SUNMemoryHelper_Sys(sunctx);
  if (!helper)
  {
    std::cout << "  SUNMemoryHelper_Sys... FAILED\n";
    return -1;
  }
  std::cout << "  SUNMemoryHelper_Sys... PASSED\n";

  test_instance(helper, SUNMEMTYPE_HOST, true);

  std::cout << "  SUNMemoryHelper_Clone... \n";
  SUNMemoryHelper helper2 = SUNMemoryHelper_Clone(helper);
  if (!helper || test_instance(helper2, SUNMEMTYPE_HOST, false))
  {
    std::cout << "  SUNMemoryHelper_Clone... FAILED\n";
    return -1;
  }
  std::cout << "  SUNMemoryHelper_Clone... PASSED\n";

  // Check destroy
  std::cout << "  SUNMemoryHelper_Destroy... \n";
  if (SUNMemoryHelper_Destroy(helper) || SUNMemoryHelper_Destroy(helper2))
  {
    std::cout << "  SUNMemoryHelper_Destroy... FAILED\n";
    return -1;
  }
  std::cout << "  SUNMemoryHelper_Destroy... PASSED\n";

  return 0;
}
