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

#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>

static std::string ReadFile(const std::string& path)
{
  std::ifstream file(path);
  std::string line;
  std::string file_contents;
  while (std::getline(file, line))
  {
    file_contents += line;
    file_contents.push_back('\n');
  }
  return file_contents;
}

static int CountLines(const std::string& s)
{
  int count = 0;
  for (char c : s)
  {
    if (c == '\n') { count++; }
  }
  return count;
}

TEST(SUNLoggerTest, EmptyFilenameDisablesWarningOutput)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_WARNING
  GTEST_SKIP() << "Warnings not enabled in this build";
#else
  const std::string warnfile = "test_sundials_logger.warn";

  (void)std::remove(warnfile.c_str());

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetWarningFilename(logger, warnfile.c_str()), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_WARNING, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(warnfile)), 1);

  ASSERT_EQ(SUNLogger_SetWarningFilename(logger, ""), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_WARNING, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(warnfile)), 1);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  (void)std::remove(warnfile.c_str());
#endif
}

TEST(SUNLoggerTest, EmptyFilenameDisablesErrorOutput)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_ERROR
  GTEST_SKIP() << "Errors not enabled in this build";
#else
  const std::string errfile = "test_sundials_logger.err";

  (void)std::remove(errfile.c_str());

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetErrorFilename(logger, errfile.c_str()), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_ERROR, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(errfile)), 1);

  ASSERT_EQ(SUNLogger_SetErrorFilename(logger, ""), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_ERROR, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(errfile)), 1);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  (void)std::remove(errfile.c_str());
#endif
}
