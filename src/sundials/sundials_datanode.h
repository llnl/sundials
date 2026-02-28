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
 * -----------------------------------------------------------------
 * SUNDataNode class definition. A SUNDataNode is a hierarchical
 * object that can hold arbitrary data in arbitrary storage locations.
 * The data may be held directly (a leaf node) or indirectly by
 * holding references to child nodes (list or object nodes). A
 * SUNDataNode maps well to a JSON node.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_DATANODE_H
#define _SUNDIALS_DATANODE_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include "sundials/sundials_memory.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef int64_t sundataindex;

enum SUNDataNodeType
{
  SUNDATANODE_LEAF,
  SUNDATANODE_LIST,
  SUNDATANODE_OBJECT
};

#ifndef SWIG
typedef enum SUNDataNodeType SUNDataNodeType;
#endif

typedef struct SUNDataNode_Ops_* SUNDataNode_Ops;
typedef struct SUNDataNode_* SUNDataNode;

struct SUNDataNode_Ops_
{
  SUNErrCode (*haschildren)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*isleaf)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*islist)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*isobject)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*addchild)(SUNDataNode, SUNDataNode child_node);
  SUNErrCode (*addnamedchild)(SUNDataNode, const char* name,
                              SUNDataNode child_node);
  SUNErrCode (*getchild)(const SUNDataNode, sundataindex index,
                         SUNDataNode* child_node);
  SUNErrCode (*getnamedchild)(const SUNDataNode, const char* name,
                              SUNDataNode* child_node);
  SUNErrCode (*removechild)(SUNDataNode, sundataindex index,
                            SUNDataNode* child_node);
  SUNErrCode (*removenamedchild)(const SUNDataNode, const char* name,
                                 SUNDataNode* child_node);
  SUNErrCode (*getdata)(const SUNDataNode, void** data, size_t* data_stride,
                        size_t* data_bytes);
  SUNErrCode (*getdatanvector)(const SUNDataNode, N_Vector v, sunrealtype* t);
  SUNErrCode (*setdata)(SUNDataNode, SUNMemoryType src_mem_type,
                        SUNMemoryType node_mem_type, void* data,
                        size_t data_stride, size_t data_bytes);
  SUNErrCode (*setdatanvector)(SUNDataNode, N_Vector v, sunrealtype t);
  SUNErrCode (*destroy)(SUNDataNode*);
};

struct SUNDataNode_
{
  SUNDataNode_Ops ops;
  SUNDataNodeType dtype;
  void* content;
  SUNContext sunctx;
};

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateLeaf(SUNDataIOMode io_mode,
                                  SUNMemoryHelper mem_helper, SUNContext sunctx,
                                  SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateList(SUNDataIOMode io_mode,
                                  sundataindex num_elements, SUNContext sunctx,
                                  SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateObject(SUNDataIOMode io_mode,
                                    sundataindex num_elements,
                                    SUNContext sunctx, SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_IsLeaf(const SUNDataNode self, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_IsList(const SUNDataNode self, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_HasChildren(const SUNDataNode self,
                                   sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_AddChild(SUNDataNode self, SUNDataNode child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode self, const char* name,
                                     SUNDataNode child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetChild(const SUNDataNode self, sundataindex index,
                                SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode self, const char* name,
                                     SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_RemoveChild(SUNDataNode self, sundataindex index,
                                   SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_RemoveNamedChild(const SUNDataNode self, const char* name,
                                        SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetData(const SUNDataNode self, void** data,
                               size_t* data_stride, size_t* data_bytes);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetDataNvector(const SUNDataNode self, N_Vector v,
                                      sunrealtype* t);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_SetData(SUNDataNode self, SUNMemoryType src_mem_type,
                               SUNMemoryType node_mem_type, void* data,
                               size_t data_stride, size_t data_bytes);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_SetDataNvector(SUNDataNode self, N_Vector v,
                                      sunrealtype t);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_Destroy(SUNDataNode* node);

#ifdef __cplusplus
}
#endif

#endif
