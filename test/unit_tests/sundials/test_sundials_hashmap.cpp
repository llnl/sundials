/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#include <iostream>
#include <limits.h>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.hpp>

#include "sundials_hashmap_impl.h"

// Helper function to free memory for value
static SUNErrCode freeKeyValue(SUNHashMapKeyValue* ptr)
{
  if (!ptr || !(*ptr)) { return SUN_SUCCESS; }
  free((*ptr)->key);
  free(*ptr);
  return SUN_SUCCESS;
}

class SUNHashMapTest : public testing::Test
{
protected:
  SUNHashMap map;

  virtual void SetUp(size_t init_capacity)
  {
    SUNHashMap_New(init_capacity, freeKeyValue, &map);
  }

  virtual void TearDown() override { SUNHashMap_Destroy(&map); }

private:
  using testing::Test::SetUp; /* silence warning from SetUp override */
};

TEST_F(SUNHashMapTest, CapacityWorks)
{
  SetUp(1);
  EXPECT_EQ(1, SUNHashMap_Capacity(map));
}

TEST_F(SUNHashMapTest, InsertAndGetWorks)
{
  SetUp(1);

  int64_t err     = 0;
  const char* key = "test_key";
  int value       = 42;

  err = SUNHashMap_Insert(map, key, &value);
  ASSERT_EQ(err, SUN_SUCCESS);

  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key, &retrieved_value);
  ASSERT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(value, *((int*)retrieved_value));
}

TEST_F(SUNHashMapTest, InsertRequiringResizeWorks)
{
  SetUp(2);

  int64_t err      = 0;
  const char* key1 = "test_key1";
  const char* key2 = "test_key2";
  const char* key3 = "test_key3";
  int value1       = 42;
  int value2       = 43;
  int value3       = 44;

  err = SUNHashMap_Insert(map, key1, &value1);
  ASSERT_EQ(err, 0);
  ASSERT_EQ(SUNHashMap_Capacity(map), 2); // Ensure no resize happened

  err = SUNHashMap_Insert(map, key2, &value2);
  ASSERT_EQ(err, 0);
  ASSERT_EQ(SUNHashMap_Capacity(map), 2); // Ensure no resize happened

  // This should trigger a resize since init_capacity is 2
  err = SUNHashMap_Insert(map, key3, &value3);
  ASSERT_EQ(err, 0);
  ASSERT_EQ(SUNHashMap_Capacity(map), 3); // Ensure resize happened

  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key1, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value1, *((int*)retrieved_value));

  err = SUNHashMap_GetValue(map, key2, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value2, *((int*)retrieved_value));

  err = SUNHashMap_GetValue(map, key3, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value3, *((int*)retrieved_value));
}

TEST_F(SUNHashMapTest, InsertDuplicateKeyFails)
{
  SetUp(1);

  int64_t err;

  // Insert same key twice (should cause error)
  const char* key = "test_key";
  int value1      = 42;
  int value2      = 100;

  err = SUNHashMap_Insert(map, key, &value1);
  ASSERT_EQ(err, 0);
  err = SUNHashMap_Insert(map, key, &value2);
  ASSERT_EQ(err, -2);
}

TEST_F(SUNHashMapTest, RemoveWorks)
{
  SetUp(2);

  int64_t err;

  // Insert a key-value pair
  const char* key = "test_key";
  int value       = 42;
  err             = SUNHashMap_Insert(map, key, &value);
  ASSERT_EQ(err, 0);

  // Remove the key
  void* removed_value;
  err = SUNHashMap_Remove(map, key, &removed_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(&value, removed_value);

  // Check if key is gone
  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key, &retrieved_value);
  ASSERT_EQ(err, -1);
}
