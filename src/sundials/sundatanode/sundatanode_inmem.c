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

#include <string.h>

#include <sundials/sundials_errors.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include "sundatanode/sundatanode_inmem.h"
#include "sundials_datanode.h"
#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"

#define GET_CONTENT(node)       ((SUNDataNode_InMemContent)(node)->content)
#define IMPL_MEMBER(node, prop) (GET_CONTENT(node)->prop)
#define BASE_MEMBER(node, prop) ((node)->prop)

static SUNDataNode sunDataNode_CreateCommon_InMem(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  SUNCheckCallNull(SUNDataNode_CreateEmpty(sunctx, &node));

  node->ops->haschildren      = SUNDataNode_HasChildren_InMem;
  node->ops->isleaf           = SUNDataNode_IsLeaf_InMem;
  node->ops->islist           = SUNDataNode_IsList_InMem;
  node->ops->isobject         = SUNDataNode_IsObject_InMem;
  node->ops->addchild         = SUNDataNode_AddChild_InMem;
  node->ops->addnamedchild    = SUNDataNode_AddNamedChild_InMem;
  node->ops->getchild         = SUNDataNode_GetChild_InMem;
  node->ops->getnamedchild    = SUNDataNode_GetNamedChild_InMem;
  node->ops->removechild      = SUNDataNode_RemoveChild_InMem;
  node->ops->removenamedchild = SUNDataNode_RemoveNamedChild_InMem;
  node->ops->getdata          = SUNDataNode_GetData_InMem;
  node->ops->getdatanvector   = SUNDataNode_GetDataNvector_InMem;
  node->ops->setdata          = SUNDataNode_SetData_InMem;
  node->ops->setdatanvector   = SUNDataNode_SetDataNvector_InMem;
  node->ops->destroy          = SUNDataNode_Destroy_InMem;

  SUNDataNode_InMemContent content =
    (SUNDataNode_InMemContent)malloc(sizeof(*content));
  SUNAssertNull(content, SUN_ERR_MEM_FAIL);

  content->parent             = NULL;
  content->mem_helper         = NULL;
  content->leaf_data          = NULL;
  content->name               = NULL;
  content->named_children     = NULL;
  content->num_named_children = 0;
  content->anon_children      = NULL;

  node->content = (void*)content;

  return node;
}

static void sunDataNode_DestroyCommon_InMem(SUNDataNode* node)
{
  if (!node || !(*node)) { return; }
  if (BASE_MEMBER(*node, content)) { free(BASE_MEMBER(*node, content)); }
  BASE_MEMBER(*node, content) = NULL;
  free(BASE_MEMBER(*node, ops));
  free(*node);
  *node = NULL;
}

static SUNErrCode sunDataNode_FreeKeyValue_InMem(SUNHashMapKeyValue* kv_ptr);
static SUNErrCode sunDataNode_FreeValue_InMem(SUNDataNode* nodeptr);

SUNErrCode SUNDataNode_CreateList_InMem(sundataindex init_size,
                                        SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNode_CreateCommon_InMem(sunctx);

  BASE_MEMBER(node, dtype) = SUNDATANODE_LIST;
  IMPL_MEMBER(node, anon_children) =
    SUNStlVector_SUNDataNode_New(init_size, sunDataNode_FreeValue_InMem);
  if (IMPL_MEMBER(node, anon_children) == NULL)
  {
    sunDataNode_DestroyCommon_InMem(&node);
    return SUN_ERR_MEM_FAIL;
  }

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateObject_InMem(sundataindex init_size,
                                          SUNContext sunctx,
                                          SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNode_CreateCommon_InMem(sunctx);

  BASE_MEMBER(node, dtype) = SUNDATANODE_OBJECT;

  SUNHashMap map;
  SUNCheckCall(SUNHashMap_New(init_size, sunDataNode_FreeKeyValue_InMem, &map));

  IMPL_MEMBER(node, named_children) = map;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf_InMem(SUNMemoryHelper mem_helper,
                                        SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNode_CreateCommon_InMem(sunctx);

  BASE_MEMBER(node, dtype)      = SUNDATANODE_LEAF;
  IMPL_MEMBER(node, mem_helper) = mem_helper;
  IMPL_MEMBER(node, leaf_data)  = NULL;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsLeaf_InMem(const SUNDataNode self,
                                    sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no = BASE_MEMBER(self, dtype) == SUNDATANODE_LEAF;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsList_InMem(const SUNDataNode self,
                                    sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no = BASE_MEMBER(self, dtype) == SUNDATANODE_LIST;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsObject_InMem(const SUNDataNode self,
                                      sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no = BASE_MEMBER(self, dtype) == SUNDATANODE_OBJECT;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_HasChildren_InMem(const SUNDataNode self,
                                         sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no =
    (IMPL_MEMBER(self, anon_children) &&
     SUNStlVector_SUNDataNode_Size(IMPL_MEMBER(self, anon_children)) != 0) ||
    IMPL_MEMBER(self, num_named_children) != 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddChild_InMem(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNAssert(BASE_MEMBER(self, dtype) == SUNDATANODE_LIST, SUN_ERR_ARG_WRONGTYPE);
  SUNStlVector_SUNDataNode_PushBack(IMPL_MEMBER(self, anon_children), child_node);
  IMPL_MEMBER(child_node, parent) = self;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddNamedChild_InMem(SUNDataNode self, const char* name,
                                           SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNAssert(BASE_MEMBER(self, dtype) == SUNDATANODE_OBJECT,
            SUN_ERR_ARG_WRONGTYPE);

  IMPL_MEMBER(child_node, name) = name;
  if (SUNHashMap_Insert(IMPL_MEMBER(self, named_children), name, child_node))
  {
    return SUN_ERR_OP_FAIL;
  }

  IMPL_MEMBER(child_node, parent) = self;
  IMPL_MEMBER(self, num_named_children)++;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetChild_InMem(const SUNDataNode self, sundataindex index,
                                      SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (!has_children) { return SUN_ERR_DATANODE_NODENOTFOUND; }

  SUNDataNode* child_node_ptr =
    SUNStlVector_SUNDataNode_At(IMPL_MEMBER(self, anon_children), index);

  if (child_node_ptr)
  {
    *child_node = *child_node_ptr;
    return SUN_SUCCESS;
  }
  else { return SUN_ERR_DATANODE_NODENOTFOUND; }
}

SUNErrCode SUNDataNode_GetNamedChild_InMem(const SUNDataNode self,
                                           const char* name,
                                           SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  *child_node = NULL;

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (has_children)
  {
    if (SUNHashMap_GetValue(IMPL_MEMBER(self, named_children), name,
                            (void**)child_node))
    {
      return SUN_ERR_DATANODE_NODENOTFOUND;
    }
    return SUN_SUCCESS;
  }
  else { return SUN_ERR_DATANODE_NODENOTFOUND; }
}

SUNErrCode SUNDataNode_RemoveChild_InMem(SUNDataNode self, sundataindex index,
                                         SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (!has_children)
  {
    *child_node = NULL;
    return SUN_SUCCESS;
  }

  SUNDataNode* child_node_ptr =
    SUNStlVector_SUNDataNode_At(IMPL_MEMBER(self, anon_children), index);
  if (child_node_ptr)
  {
    *child_node = *child_node_ptr;
    if (*child_node)
    {
      IMPL_MEMBER(*child_node, parent) = NULL;
      SUNStlVector_SUNDataNode_Erase(IMPL_MEMBER(self, anon_children), index);
    }
    else { return SUN_ERR_DATANODE_NODENOTFOUND; }
  }
  else { return SUN_ERR_DATANODE_NODENOTFOUND; }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_RemoveNamedChild_InMem(const SUNDataNode self,
                                              const char* name,
                                              SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  *child_node = NULL;

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (has_children)
  {
    if (SUNHashMap_Remove(IMPL_MEMBER(self, named_children), name,
                          (void**)child_node))
    {
      *child_node = NULL;
      return SUN_ERR_DATANODE_NODENOTFOUND;
    }
    IMPL_MEMBER(*child_node, parent) = NULL;
    IMPL_MEMBER(self, num_named_children)--;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetData_InMem(const SUNDataNode self, void** data,
                                     size_t* data_stride, size_t* data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNMemory leaf_data = (SUNMemory)IMPL_MEMBER(self, leaf_data);

  *data_stride = leaf_data->stride;
  *data_bytes  = leaf_data->bytes;
  *data        = leaf_data->ptr;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetDataNvector_InMem(const SUNDataNode self, N_Vector v,
                                            sunrealtype* t)
{
  SUNFunctionBegin(self->sunctx);

  /* Use the default queue for the memory helper */
  void* queue = NULL;

  SUNMemory leaf_data = (SUNMemory)IMPL_MEMBER(self, leaf_data);

  SUNMemoryType leaf_mem_type = leaf_data->type;

  sunindextype buffer_size = 0;
  SUNCheckCall(N_VBufSize(v, &buffer_size));
  SUNAssert((buffer_size + sizeof(sunrealtype)) == leaf_data->bytes,
            SUN_ERR_ARG_INCOMPATIBLE);

  if (leaf_mem_type != SUNMEMTYPE_HOST)
  {
    /* BufUnpack assumes the data is on the host. So if the leaf has it elsewhere,
       we need to move it to the host first. */
    SUNMemory leaf_host_data = NULL;
    SUNCheckCall(SUNMemoryHelper_Alloc(IMPL_MEMBER(self, mem_helper),
                                       &leaf_host_data, leaf_data->bytes,
                                       SUNMEMTYPE_HOST, queue));

    SUNCheckCall(SUNMemoryHelper_Copy(IMPL_MEMBER(self, mem_helper),
                                      leaf_host_data, leaf_data, buffer_size,
                                      queue));

    sunrealtype* data_ptr = leaf_host_data->ptr;
    *t                    = data_ptr[0];
    SUNCheckCall(N_VBufUnpack(v, &data_ptr[1]));

    SUNCheckCall(SUNMemoryHelper_Dealloc(IMPL_MEMBER(self, mem_helper),
                                         leaf_host_data, queue));
  }
  else
  {
    sunrealtype* data_ptr = leaf_data->ptr;
    *t                    = data_ptr[0];
    SUNCheckCall(N_VBufUnpack(v, &data_ptr[1]));
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetData_InMem(SUNDataNode self, SUNMemoryType src_mem_type,
                                     SUNMemoryType node_mem_type, void* data,
                                     size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  /* Use the default queue for the memory helper */
  void* queue = NULL;

  SUNAssert(BASE_MEMBER(self, dtype) == SUNDATANODE_LEAF, SUN_ERR_ARG_WRONGTYPE);

  SUNMemory data_mem_src = SUNMemoryHelper_Wrap(IMPL_MEMBER(self, mem_helper),
                                                data, src_mem_type);
  SUNCheckLastErr();

  SUNMemory data_mem_dst = NULL;
  SUNCheckCall(SUNMemoryHelper_AllocStrided(IMPL_MEMBER(self, mem_helper),
                                            &data_mem_dst, data_bytes,
                                            data_stride, node_mem_type, queue));

  SUNCheckCall(SUNMemoryHelper_Copy(IMPL_MEMBER(self, mem_helper), data_mem_dst,
                                    data_mem_src, data_bytes, queue));

  SUNMemoryHelper_Dealloc(IMPL_MEMBER(self, mem_helper), data_mem_src, queue);

  IMPL_MEMBER(self, leaf_data) = data_mem_dst;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetDataNvector_InMem(SUNDataNode self, N_Vector v,
                                            sunrealtype t)
{
  SUNFunctionBegin(self->sunctx);

  /* Use the default queue for the memory helper */
  void* queue = NULL;

  SUNMemoryType leaf_mem_type = SUNMEMTYPE_HOST;

  sunindextype buffer_size = 0;
  SUNCheckCall(N_VBufSize(v, &buffer_size));

  /* We allocate 1 extra sunrealtype for storing t */
  SUNMemory leaf_data = NULL;
  SUNCheckCall(
    SUNMemoryHelper_AllocStrided(IMPL_MEMBER(self, mem_helper), &leaf_data,
                                 buffer_size + sizeof(sunrealtype),
                                 sizeof(sunrealtype), leaf_mem_type, queue));

  /* BufPack will handle any necessary copies from the device and will fill data_ptr on the host */
  sunrealtype* data_ptr = leaf_data->ptr;
  data_ptr[0]           = t;
  SUNCheckCall(N_VBufPack(v, &data_ptr[1]));

  IMPL_MEMBER(self, leaf_data) = leaf_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_Destroy_InMem(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  /* Use the default queue for the memory helper */
  void* queue = NULL;

  if (BASE_MEMBER(*node, dtype) == SUNDATANODE_OBJECT)
  {
    SUNHashMap map = IMPL_MEMBER(*node, named_children);
    SUNHashMap_Destroy(&map);
  }
  else if (BASE_MEMBER(*node, dtype) == SUNDATANODE_LIST)
  {
    SUNStlVector_SUNDataNode_Destroy(&IMPL_MEMBER(*node, anon_children));
  }
  else if (BASE_MEMBER(*node, dtype) == SUNDATANODE_LEAF)
  {
    if (IMPL_MEMBER(*node, leaf_data))
    {
      SUNCheckCall(SUNMemoryHelper_Dealloc(IMPL_MEMBER(*node, mem_helper),
                                           IMPL_MEMBER(*node, leaf_data), queue));
    }
  }

  sunDataNode_DestroyCommon_InMem(node);
  *node = NULL;

  return SUN_SUCCESS;
}

/* This function is the callback provided to the child hashmap as the destroy function. */
static SUNErrCode sunDataNode_FreeKeyValue_InMem(SUNHashMapKeyValue* kv_ptr)
{
  if (!kv_ptr || !(*kv_ptr)) { return SUN_SUCCESS; }
  SUNDataNode node = (SUNDataNode)((*kv_ptr)->value);
  free((*kv_ptr)->key);
  free(*kv_ptr);
  SUNDataNode_Destroy_InMem(&node);
  return SUN_SUCCESS;
}

/* This function is the callback provided to the child stlvector as the destroy function. */
static SUNErrCode sunDataNode_FreeValue_InMem(SUNDataNode* nodeptr)
{
  SUNDataNode_Destroy_InMem(nodeptr);
  return SUN_SUCCESS;
}
