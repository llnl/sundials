/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * SUNAdjointCheckpointScheme private class definition.
 * ----------------------------------------------------------------*/

#ifndef SUNDIALS_ADJOINTCHECKPOINTSCHEME_IMPL_H_
#define SUNDIALS_ADJOINTCHECKPOINTSCHEME_IMPL_H_

#include <sundials/sundials_adjointcheckpointscheme.h>

typedef struct SUNAdjointCheckpointScheme_Ops_* SUNAdjointCheckpointScheme_Ops;

struct SUNAdjointCheckpointScheme_Ops_
{
  SUNAdjointCheckpointSchemeNeedsSavingFn needssaving;
  SUNAdjointCheckpointSchemeInsertVectorFn insertvector;
  SUNAdjointCheckpointSchemeLoadVectorFn loadvector;
  SUNAdjointCheckpointSchemeDestroyFn destroy;
  SUNAdjointCheckpointSchemeEnableDenseFn enableDense;
};

struct SUNAdjointCheckpointScheme_
{
  SUNAdjointCheckpointScheme_Ops ops;
  void* content;
  SUNContext sunctx;
};

#endif /* SUNDIALS_ADJOINTCHECKPOINTSCHEME_IMPL_H_ */
