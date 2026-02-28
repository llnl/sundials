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

#include <arkode/arkode.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"

static const std::string errfile{"test_arkode_error_handling.err"};

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  return 0;
}

class ARKodeErrConditionTest : public testing::Test
{
protected:
  ARKodeErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    v          = N_VNew_Serial(1, sunctx);
    arkode_mem = ARKStepCreate(f, f, sunrealtype{0.0}, v, sunctx);
  }

  ~ARKodeErrConditionTest()
  {
    N_VDestroy(v);
    ARKodeFree(&arkode_mem);
  }

  void* arkode_mem;
  N_Vector v;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(ARKodeErrConditionTest, WarningIsPrinted)
{
  SUNErrCode err = SUNLogger_SetWarningFilename(logger, errfile.c_str());
  ASSERT_EQ(err, SUN_SUCCESS);
  ARKodeMemRec* ark_mem = (ARKodeMemRec*)arkode_mem;
  arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__, "test");
  err = SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING);
  ASSERT_EQ(err, SUN_SUCCESS);
  std::string output = dumpstderr(sunctx, errfile);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[WARNING]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr("test")));
#else
  EXPECT_EQ(output, "");
#endif
}

TEST_F(ARKodeErrConditionTest, ErrorIsPrinted)
{
  SUNErrCode err = SUNLogger_SetErrorFilename(logger, errfile.c_str());
  ASSERT_EQ(err, SUN_SUCCESS);
  // negative reltol is illegal
  int ierr = ARKodeSStolerances(arkode_mem, /* reltol= */ -1e-4,
                                /* abstol= */ 1e-4);
  ASSERT_NE(ierr, ARK_SUCCESS);
  std::string output = dumpstderr(sunctx, errfile);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_ERROR
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSG_ARK_BAD_RELTOL)));
#else
  EXPECT_EQ(output, "");
#endif
}
