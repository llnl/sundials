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

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.h>
#include <sunmemory/sunmemory_system.h>
#include "sundials_datanode.h"

#include "sundatanode/sundatanode_inmem.h"
#include "sundials/sundials_memory.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

#define GET_IMPL(node)       ((SUNDataNode_InMemContent)(node)->content)
#define GET_PROP(node, prop) (GET_IMPL(node)->prop)

static int get_leaf_as_int(SUNDataNode node)
{
  SUNMemory mem = (SUNMemory)GET_PROP(node, leaf_data);
  return *((int*)mem->ptr);
}

class SUNDataNodeTest : public testing::Test
{
protected:
  SUNDataNodeTest()
  {
    SUNContext_Create(SUN_COMM_NULL, &sunctx);
    mem_helper = SUNMemoryHelper_Sys(sunctx);
  }

  ~SUNDataNodeTest()
  {
    SUNContext_Free(&sunctx);
    SUNMemoryHelper_Destroy_Sys(mem_helper);
  }

  SUNContext sunctx;
  SUNMemoryHelper mem_helper;
};

TEST_F(SUNDataNodeTest, CreateEmptyWorks)
{
  SUNErrCode err;
  SUNDataNode node;

  err = SUNDataNode_CreateEmpty(sunctx, &node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, CreateLeafWorks)
{
  SUNErrCode err;
  SUNDataNode node;
  int integer_value = 5;

  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx, &node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(integer_value, get_leaf_as_int(node));

  err = SUNDataNode_Destroy(&node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, CreateLeafWorksWhenEmpty)
{
  SUNErrCode err;
  SUNDataNode node;

  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx, &node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(NULL, GET_PROP(node, leaf_data));

  err = SUNDataNode_Destroy(&node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, AddChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(root_node, GET_PROP(child_node, parent));

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, AddNamedChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateObject(SUNDATAIOMODE_INMEM, num_elem, sunctx,
                                 &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddNamedChild(root_node, "int_value", child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(1, GET_PROP(root_node, num_named_children));
  EXPECT_EQ(root_node, GET_PROP(child_node, parent));

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, HasChildrenWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunbooleantype yes_or_no = SUNTRUE;
  err                      = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_FALSE(yes_or_no);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_TRUE(yes_or_no);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, RemoveChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_RemoveChild(root_node, 0, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunbooleantype yes_or_no = SUNTRUE;
  err                      = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_FALSE(yes_or_no);

  EXPECT_FALSE(GET_IMPL(child_node)->parent);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, RemoveSameChildTwiceWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node, child_node_out;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_RemoveChild(root_node, 0, &child_node_out);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(child_node_out, child_node);

  err = SUNDataNode_RemoveChild(root_node, 0, &child_node_out);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, RemoveChildWorksWhenEmpty)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node, child_node_out;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_RemoveChild(root_node, 0, &child_node_out);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, GetChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_GetChild(root_node, 0, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(integer_value, get_leaf_as_int(child_node));

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, GetNamedChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateObject(SUNDATAIOMODE_INMEM, num_elem, sunctx,
                                 &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddNamedChild(root_node, "int_value", child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_GetNamedChild(root_node, "int_value", &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  size_t stride, bytes;
  err = SUNDataNode_GetData(child_node, &raw_value, &stride, &bytes);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(integer_value, *static_cast<int*>(raw_value));

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, RemoveNamedChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node, child_node_out;
  unsigned int num_elem = 5;

  err = SUNDataNode_CreateObject(SUNDATAIOMODE_INMEM, num_elem, sunctx,
                                 &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(child_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddNamedChild(root_node, "int_value", child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_RemoveNamedChild(root_node, "int_value", &child_node_out);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(integer_value, get_leaf_as_int(child_node_out));

  sunbooleantype yes_or_no = SUNTRUE;
  err                      = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_FALSE(yes_or_no);
  EXPECT_FALSE(GET_IMPL(child_node)->parent);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, GetDataWorksWhenLeaf)
{
  SUNErrCode err;
  SUNDataNode root_node;

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_SetData(root_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)(&integer_value), sizeof(integer_value),
                            sizeof(integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  size_t stride, bytes;
  err = SUNDataNode_GetData(root_node, &raw_value, &stride, &bytes);
  EXPECT_EQ(err, SUN_SUCCESS);

  int value_we_got = *((int*)raw_value);
  EXPECT_EQ(integer_value, value_we_got);
  EXPECT_EQ(stride, sizeof(integer_value));
  EXPECT_EQ(bytes, sizeof(integer_value));

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, SetDataWorksWhenLeaf)
{
  SUNErrCode err;
  SUNDataNode root_node;

  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int new_integer_value = 3;
  err = SUNDataNode_SetData(root_node, SUNMEMTYPE_HOST, SUNMEMTYPE_HOST,
                            (void*)&new_integer_value, sizeof(new_integer_value),
                            sizeof(new_integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  size_t stride, bytes;
  err = SUNDataNode_GetData(root_node, &raw_value, &stride, &bytes);
  EXPECT_EQ(err, SUN_SUCCESS);

  int value_we_got = *((int*)raw_value);
  EXPECT_EQ(new_integer_value, value_we_got);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, SetAndGetDataNvectorWhenLeaf)
{
  SUNErrCode err;
  SUNDataNode root_node;
  sunrealtype real_value = SUN_RCONST(3.0);
  N_Vector v             = N_VNew_Serial(2, sunctx);
  N_Vector vec_we_got    = N_VClone(v);

  N_VConst(real_value, v);
  N_VConst(SUN_RCONST(0.0), vec_we_got);

  err = SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, mem_helper, sunctx,
                               &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunrealtype t = SUN_RCONST(1.0);
  err           = SUNDataNode_SetDataNvector(root_node, v, t);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunrealtype tout = SUN_RCONST(0.0);
  err              = SUNDataNode_GetDataNvector(root_node, vec_we_got, &tout);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(t, tout);
  EXPECT_EQ(N_VGetArrayPointer(v)[0], N_VGetArrayPointer(vec_we_got)[0]);
  EXPECT_EQ(N_VGetArrayPointer(v)[1], N_VGetArrayPointer(vec_we_got)[1]);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  N_VDestroy(v);
  N_VDestroy(vec_we_got);
}
