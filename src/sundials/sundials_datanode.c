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
 * -----------------------------------------------------------------*
 * The SUNDataNode class is a hierarchical object which could be used
 * to build something like a JSON tree. Nodes can be lists, objects,
 * or leaves. Using the JSON analogy:
 *
 *   "i_am_object": {
 *     "i_am_list": [
 *       "i_am_object": {...
 *       },
 *       "i_am_leaf"
 *     ],
 *     "i_am_leaf": "value"
 *   }
 *
 * Object nodes hold named nodes (children), while list nodes hold
 * anonymous nodes (children). Leaf nodes do not have children, they
 * have values. The SUNDataNode can be used to build all sorts of
 * useful things, but we primarily use it as the backbone for
 * checkpointing states in adjoint sensitivity analysis.
 * -----------------------------------------------------------------*/

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundatanode/sundatanode_inmem.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_memory.h"
#include "sundials_datanode.h"

/**
 * :param sunctx: The SUNContext.
 * :param node_out: Pointer to the output SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode self;
  self = (SUNDataNode)malloc(sizeof(*self));
  SUNAssert(self, SUN_ERR_MEM_FAIL);

  SUNDataNode_Ops ops;
  ops = (SUNDataNode_Ops)malloc(sizeof(*ops));
  SUNAssert(self, SUN_ERR_MEM_FAIL);

  ops->haschildren = NULL;
  ops->isleaf      = NULL;
  ops->islist      = NULL;
  ops->isobject    = NULL;
  ops->addchild    = NULL;
  ops->getchild    = NULL;
  ops->removechild = NULL;
  ops->getdata     = NULL;
  ops->setdata     = NULL;
  ops->destroy     = NULL;

  self->dtype   = 0;
  self->ops     = ops;
  self->content = NULL;
  self->sunctx  = sunctx;

  *node_out = self;
  return SUN_SUCCESS;
}

/**
 * :param io_mode: The I/O mode used for storing the data.
 * :param mem_helper: The memory helper.
 * :param sunctx: The SUNContext.
 * :param node_out: Pointer to the output SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_CreateLeaf(SUNDataIOMode io_mode,
                                  SUNMemoryHelper mem_helper, SUNContext sunctx,
                                  SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  SUNErrCode err = SUN_SUCCESS;
  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    err = SUNDataNode_CreateLeaf_InMem(mem_helper, sunctx, node_out);
    break;
  default: err = SUN_ERR_ARG_OUTOFRANGE;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  SUNCheck(err == SUN_SUCCESS, err);

  return err;
}

/**
 * :param io_mode: The I/O mode used for storing the data.
 * :param num_elements: The number of elements in the list.
 * :param sunctx: The SUNContext.
 * :param node_out: Pointer to the output SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_CreateList(SUNDataIOMode io_mode,
                                  sundataindex num_elements, SUNContext sunctx,
                                  SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  SUNErrCode err = SUN_SUCCESS;
  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    err = SUNDataNode_CreateList_InMem(num_elements, sunctx, node_out);
    break;
  default: err = SUN_ERR_ARG_OUTOFRANGE;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  SUNCheck(err == SUN_SUCCESS, err);

  return err;
}

/**
 * :param io_mode: The I/O mode used for storing the data.
 * :param num_elements: The number of elements in the object.
 * :param sunctx: The SUNContext.
 * :param node_out: Pointer to the output SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_CreateObject(SUNDataIOMode io_mode,
                                    sundataindex num_elements,
                                    SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  SUNErrCode err = SUN_SUCCESS;
  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    err = SUNDataNode_CreateObject_InMem(num_elements, sunctx, node_out);
    break;
  default: err = SUN_ERR_ARG_OUTOFRANGE;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  SUNCheck(err == SUN_SUCCESS, err);

  return err;
}

/**
 * :param self: The SUNDataNode.
 * :param yes_or_no: Pointer to the output boolean result.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_IsLeaf(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->isleaf)
  {
    SUNErrCode err = self->ops->isleaf(self, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param yes_or_no: Pointer to the output boolean result.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_IsList(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->islist)
  {
    SUNErrCode err = self->ops->islist(self, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param yes_or_no: Pointer to the output boolean result.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_HasChildren(const SUNDataNode self,
                                   sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->haschildren)
  {
    SUNErrCode err = self->ops->haschildren(self, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param child_node: The child SUNDataNode to add.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_AddChild(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->addchild)
  {
    SUNErrCode err = self->ops->addchild(self, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param name: The name of the child.
 * :param child_node: The child SUNDataNode to add.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode self, const char* name,
                                     SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->addnamedchild)
  {
    SUNErrCode err = self->ops->addnamedchild(self, name, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param index: The index of the child.
 * :param child_node: Pointer to the output child SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_GetChild(const SUNDataNode self, sundataindex index,
                                SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getchild)
  {
    SUNErrCode err = self->ops->getchild(self, index, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param name: The name of the child.
 * :param child_node: Pointer to the output child SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode self, const char* name,
                                     SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getnamedchild)
  {
    SUNErrCode err = self->ops->getnamedchild(self, name, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param name: The name of the child.
 * :param child_node: Pointer to the output child SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_RemoveNamedChild(const SUNDataNode self, const char* name,
                                        SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->removenamedchild)
  {
    SUNErrCode err = self->ops->removenamedchild(self, name, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param index: The index of the child.
 * :param child_node: Pointer to the output child SUNDataNode.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_RemoveChild(SUNDataNode self, sundataindex index,
                                   SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->removechild)
  {
    SUNErrCode err = self->ops->removechild(self, index, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param data: Pointer to the output data.
 * :param data_stride: Pointer to the output data stride.
 * :param data_bytes: Pointer to the output data bytes.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_GetData(const SUNDataNode self, void** data,
                               size_t* data_stride, size_t* data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getdata)
  {
    SUNErrCode err = self->ops->getdata(self, data, data_stride, data_bytes);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param v: The output state N_Vector.
 * :param t: On output, the time associated with the output state vector.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_GetDataNvector(SUNDataNode self, N_Vector v, sunrealtype* t)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getdatanvector)
  {
    SUNErrCode err = self->ops->getdatanvector(self, v, t);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param src_mem_type: The source memory type.
 * :param node_mem_type: The node memory type.
 * :param data: The data to set.
 * :param data_stride: The data stride.
 * :param data_bytes: The data bytes.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_SetData(SUNDataNode self, SUNMemoryType src_mem_type,
                               SUNMemoryType node_mem_type, void* data,
                               size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->setdata)
  {
    SUNErrCode err = self->ops->setdata(self, src_mem_type, node_mem_type, data,
                                        data_stride, data_bytes);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param self: The SUNDataNode.
 * :param v: The state N_Vector.
 * :param t: The time associated with the state vector.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_SetDataNvector(SUNDataNode self, N_Vector v, sunrealtype t)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->setdatanvector)
  {
    SUNErrCode err = self->ops->setdatanvector(self, v, t);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

/**
 * :param node: Pointer to the SUNDataNode to destroy.
 * :return: SUNErrCode indicating success or failure.
 */
SUNErrCode SUNDataNode_Destroy(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if ((*node)->ops->destroy)
  {
    SUNErrCode err = (*node)->ops->destroy(node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  free((*node)->ops);
  free(*node);
  *node = NULL;

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_SUCCESS;
}
