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

#include <cvode/cvode.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "cvode/cvode_impl.h"
#include "sundials/sundials_context.hpp"

static const std::string errfile{"test_error_handling.err"};

class CVodeErrConditionTest : public testing::Test
{
protected:
  CVodeErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
  }

  ~CVodeErrConditionTest() { CVodeFree(&cvode_mem); }

  void* cvode_mem;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(CVodeErrConditionTest, WarningIsPrinted)
{
  SUNErrCode err = SUNLogger_SetWarningFilename(logger, errfile.c_str());
  ASSERT_EQ(err, SUN_SUCCESS);
  CVodeMemRec* ark_mem = (CVodeMemRec*)cvode_mem;
  cvProcessError(ark_mem, CV_WARNING, __LINE__, __func__, __FILE__, "test");
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

TEST_F(CVodeErrConditionTest, ErrorIsPrinted)
{
  SUNErrCode err = SUNLogger_SetErrorFilename(logger, errfile.c_str());
  ASSERT_EQ(err, SUN_SUCCESS);
  // attempting to call CVodeSStolerances before CVodeInit is illegal
  int ierr = CVodeSStolerances(cvode_mem, 1e-4, 1e-4);
  ASSERT_NE(ierr, CV_SUCCESS);
  std::string output = dumpstderr(sunctx, errfile);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_ERROR
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSGCV_NO_MALLOC)));
#else
  EXPECT_EQ(output, "");
#endif
}
