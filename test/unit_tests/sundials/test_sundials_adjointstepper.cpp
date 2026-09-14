/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, Lawrence Livermore National Security,
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

#include <gtest/gtest.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_stepper.h>

#include "sundials_adjointstepper_impl.h"

struct ReInitData
{
  int calls          = 0;
  sunrealtype t0     = SUN_RCONST(0.0);
  N_Vector y0        = nullptr;
  SUNErrCode retcode = SUN_SUCCESS;
};

static SUNErrCode reinit(SUNStepper stepper, sunrealtype t0, N_Vector y0)
{
  void* content  = nullptr;
  SUNErrCode err = SUNStepper_GetContent(stepper, &content);
  if (err != SUN_SUCCESS) { return err; }

  auto data = static_cast<ReInitData*>(content);
  data->calls++;
  data->t0 = t0;
  data->y0 = y0;
  return data->retcode;
}

TEST(SUNAdjointStepper, ReInit)
{
  SUNContext sunctx = nullptr;
  ASSERT_EQ(SUNContext_Create(SUN_COMM_NULL, &sunctx), SUN_SUCCESS);

  N_Vector sf = N_VNew_Serial(1, sunctx);
  ASSERT_NE(sf, nullptr);

  SUNStepper fwd_stepper = nullptr;
  SUNStepper adj_stepper = nullptr;
  ASSERT_EQ(SUNStepper_Create(sunctx, &fwd_stepper), SUN_SUCCESS);
  ASSERT_EQ(SUNStepper_Create(sunctx, &adj_stepper), SUN_SUCCESS);

  ReInitData fwd_data;
  ReInitData adj_data;
  ASSERT_EQ(SUNStepper_SetContent(fwd_stepper, &fwd_data), SUN_SUCCESS);
  ASSERT_EQ(SUNStepper_SetContent(adj_stepper, &adj_data), SUN_SUCCESS);
  ASSERT_EQ(SUNStepper_SetReInitFn(fwd_stepper, reinit), SUN_SUCCESS);
  ASSERT_EQ(SUNStepper_SetReInitFn(adj_stepper, reinit), SUN_SUCCESS);

  SUNAdjointCheckpointScheme checkpoint_scheme = nullptr;
  ASSERT_EQ(SUNAdjointCheckpointScheme_NewEmpty(sunctx, &checkpoint_scheme),
            SUN_SUCCESS);

  SUNAdjointStepper self = nullptr;
  ASSERT_EQ(SUNAdjointStepper_Create(fwd_stepper, SUNFALSE, adj_stepper,
                                     SUNFALSE, 2, SUN_RCONST(1.0), sf,
                                     checkpoint_scheme, sunctx, &self),
            SUN_SUCCESS);

  self->nrecompute = 3;
  ASSERT_EQ(SUNAdjointStepper_ReInit(self, SUN_RCONST(2.0), sf, 7), SUN_SUCCESS);
  EXPECT_EQ(fwd_data.calls, 0);
  EXPECT_EQ(adj_data.calls, 1);
  EXPECT_EQ(adj_data.t0, SUN_RCONST(2.0));
  EXPECT_EQ(adj_data.y0, sf);
  EXPECT_EQ(self->tf, SUN_RCONST(2.0));
  EXPECT_EQ(self->final_step_idx, 7);
  EXPECT_EQ(self->nrecompute, 0);

  self->nrecompute = 4;
  adj_data.retcode = SUN_ERR_OP_FAIL;
  EXPECT_EQ(SUNAdjointStepper_ReInit(self, SUN_RCONST(3.0), sf, 9),
            SUN_ERR_OP_FAIL);
  EXPECT_EQ(self->tf, SUN_RCONST(2.0));
  EXPECT_EQ(self->final_step_idx, 7);
  EXPECT_EQ(self->nrecompute, 4);

  EXPECT_EQ(SUNAdjointStepper_Destroy(&self), SUN_SUCCESS);
  EXPECT_EQ(SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme), SUN_SUCCESS);
  EXPECT_EQ(SUNStepper_Destroy(&adj_stepper), SUN_SUCCESS);
  EXPECT_EQ(SUNStepper_Destroy(&fwd_stepper), SUN_SUCCESS);
  N_VDestroy(sf);
  EXPECT_EQ(SUNContext_Free(&sunctx), SUN_SUCCESS);
}
