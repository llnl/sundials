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

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_types.h>

#include "sundials_datanode.h"
#include "sundials_hashmap_impl.h"

#ifndef _SUNDATANODE_INMEM_H
#define _SUNDATANODE_INMEM_H

#ifdef __cplusplus
extern "C" {
#endif

#define TTYPE SUNDataNode
#include "stl/sunstl_vector.h"
#undef TTYPE

typedef struct SUNDataNode_InMemContent_* SUNDataNode_InMemContent;

struct SUNDataNode_InMemContent_
{
  // Reference to the parent node of this node.
  SUNDataNode parent;

  // Node can only be an object, leaf, or list. It cannot be more than one of these at a time.

  // Properties for Leaf nodes (nodes that store data)
  SUNMemoryHelper mem_helper;
  SUNMemory leaf_data;

  // Properties for Object nodes (nodes that are a collection of named nodes)
  const char* name;
  SUNHashMap named_children;
  sundataindex num_named_children;

  // Properties for a List node (nodes that are a collection of anonymous nodes)
  SUNStlVector_SUNDataNode anon_children;
};

SUNErrCode SUNDataNode_CreateList_InMem(sundataindex init_size,
                                        SUNContext sunctx, SUNDataNode* node_out);

SUNErrCode SUNDataNode_CreateObject_InMem(sundataindex init_size,
                                          SUNContext sunctx,
                                          SUNDataNode* node_out);

SUNErrCode SUNDataNode_CreateLeaf_InMem(SUNMemoryHelper mem_helper,
                                        SUNContext sunctx, SUNDataNode* node_out);

SUNErrCode SUNDataNode_IsLeaf_InMem(const SUNDataNode self,
                                    sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_IsList_InMem(const SUNDataNode self,
                                    sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_IsObject_InMem(const SUNDataNode self,
                                      sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_HasChildren_InMem(const SUNDataNode self,
                                         sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_AddChild_InMem(SUNDataNode parent_node,
                                      SUNDataNode child_node);

SUNErrCode SUNDataNode_AddNamedChild_InMem(SUNDataNode parent_node,
                                           const char* name,
                                           SUNDataNode child_node);

SUNErrCode SUNDataNode_GetChild_InMem(const SUNDataNode self, sundataindex index,
                                      SUNDataNode* child_node);

SUNErrCode SUNDataNode_GetNamedChild_InMem(const SUNDataNode self,
                                           const char* name,
                                           SUNDataNode* child_node);

SUNErrCode SUNDataNode_RemoveChild_InMem(SUNDataNode self, sundataindex index,
                                         SUNDataNode* child_node);

SUNErrCode SUNDataNode_RemoveNamedChild_InMem(SUNDataNode self, const char* name,
                                              SUNDataNode* child_node);

SUNErrCode SUNDataNode_GetData_InMem(const SUNDataNode self, void** data,
                                     size_t* data_stride, size_t* data_bytes);

SUNErrCode SUNDataNode_GetDataNvector_InMem(const SUNDataNode self, N_Vector v,
                                            sunrealtype* t);

SUNErrCode SUNDataNode_SetData_InMem(SUNDataNode self, SUNMemoryType src_mem_type,
                                     SUNMemoryType node_mem_type, void* data,
                                     size_t data_stride, size_t data_bytes);

SUNErrCode SUNDataNode_SetDataNvector_InMem(SUNDataNode self, N_Vector v,
                                            sunrealtype t);

SUNErrCode SUNDataNode_Destroy_InMem(SUNDataNode* node);

#ifdef __cplusplus
}
#endif

#endif // _SUNDATANODE_INMEM_H
