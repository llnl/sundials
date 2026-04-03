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

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

#include <sundials/priv/sundials_logger_macros.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>

// Add [[maybe_unused]] because ReadFile and CountLines are "unused" when the
// logging level is set to 0

[[maybe_unused]] static std::string ReadFile(const std::string& path)
{
  std::ifstream file(path);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

[[maybe_unused]] static int CountLines(const std::string& s)
{
  return std::count(s.begin(), s.end(), '\n');
}

TEST(SUNLoggerTest, DefaultErrorLogger)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_ERROR
  GTEST_SKIP() << "Errors not enabled in this build";
#else
  testing::internal::CaptureStderr();

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_ERROR, "scope", "label",
                               "test"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);

  EXPECT_EQ(CountLines(testing::internal::GetCapturedStderr()), 1);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  ASSERT_EQ(logger, nullptr);
#endif
}

TEST(SUNLoggerTest, DefaultWarningLogger)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_WARNING
  GTEST_SKIP() << "Warnings not enabled in this build";
#else
  testing::internal::CaptureStdout();

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_WARNING, "scope", "label",
                               "test"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);

  EXPECT_EQ(CountLines(testing::internal::GetCapturedStdout()), 1);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  ASSERT_EQ(logger, nullptr);
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

TEST(SUNLoggerTest, EmptyFilenameDisablesInfoOutput)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_INFO
  GTEST_SKIP() << "Info logging not enabled in this build";
#else
  const std::string infofile = "test_sundials_logger.info";

  (void)std::remove(infofile.c_str());

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetInfoFilename(logger, infofile.c_str()), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_INFO, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_INFO), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(infofile)), 1);

  ASSERT_EQ(SUNLogger_SetInfoFilename(logger, ""), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_INFO, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_INFO), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(infofile)), 1);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  (void)std::remove(infofile.c_str());
#endif
}

TEST(SUNLoggerTest, EmptyFilenameDisablesDebugOutput)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_DEBUG
  GTEST_SKIP() << "Debug logging not enabled in this build";
#else
  const std::string debugfile = "test_sundials_logger.debug";

  (void)std::remove(debugfile.c_str());

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetDebugFilename(logger, debugfile.c_str()), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_DEBUG), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(debugfile)), 1);

  ASSERT_EQ(SUNLogger_SetDebugFilename(logger, ""), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_DEBUG), SUN_SUCCESS);

  EXPECT_EQ(CountLines(ReadFile(debugfile)), 1);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  (void)std::remove(debugfile.c_str());
#endif
}

TEST(SUNLoggerTest, SetFileNullLogger)
{
  EXPECT_EQ(SUNLogger_SetErrorFile(NULL, NULL), SUN_ERR_ARG_CORRUPT);
  EXPECT_EQ(SUNLogger_SetWarningFile(NULL, NULL), SUN_ERR_ARG_CORRUPT);
  EXPECT_EQ(SUNLogger_SetInfoFile(NULL, NULL), SUN_ERR_ARG_CORRUPT);
  EXPECT_EQ(SUNLogger_SetDebugFile(NULL, NULL), SUN_ERR_ARG_CORRUPT);
}

TEST(SUNLoggerTest, SetErrorFileAndDisableWithNull)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_ERROR
  GTEST_SKIP() << "Errors not enabled in this build";
#else
  const std::string errfile = "test_sundials_logger.error_fp.err";
  (void)std::remove(errfile.c_str());

  FILE* fp = std::fopen(errfile.c_str(), "w+");
  ASSERT_NE(fp, nullptr);

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetErrorFile(logger, fp), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_ERROR, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetErrorFile(logger, NULL), SUN_SUCCESS);
  testing::internal::CaptureStderr();
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_ERROR, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);
  EXPECT_EQ(CountLines(testing::internal::GetCapturedStderr()), 0);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  ASSERT_EQ(logger, nullptr);

  ASSERT_EQ(std::fclose(fp), 0);
  EXPECT_EQ(CountLines(ReadFile(errfile)), 1);
  (void)std::remove(errfile.c_str());
#endif
}

TEST(SUNLoggerTest, SetWarningFileAndDisableWithNull)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_WARNING
  GTEST_SKIP() << "Warnings not enabled in this build";
#else
  const std::string warnfile = "test_sundials_logger.warning_fp.warn";
  (void)std::remove(warnfile.c_str());

  FILE* fp = std::fopen(warnfile.c_str(), "w+");
  ASSERT_NE(fp, nullptr);

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetWarningFile(logger, fp), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_WARNING, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetWarningFile(logger, NULL), SUN_SUCCESS);
  testing::internal::CaptureStdout();
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_WARNING, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING), SUN_SUCCESS);
  EXPECT_EQ(CountLines(testing::internal::GetCapturedStdout()), 0);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  ASSERT_EQ(logger, nullptr);

  ASSERT_EQ(std::fclose(fp), 0);
  EXPECT_EQ(CountLines(ReadFile(warnfile)), 1);
  (void)std::remove(warnfile.c_str());
#endif
}

TEST(SUNLoggerTest, SetInfoFileAndDisableWithNull)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_INFO
  GTEST_SKIP() << "Info logging not enabled in this build";
#else
  const std::string infofile = "test_sundials_logger.info_fp.info";
  (void)std::remove(infofile.c_str());

  FILE* fp = std::fopen(infofile.c_str(), "w+");
  ASSERT_NE(fp, nullptr);

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetInfoFile(logger, fp), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_INFO, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_INFO), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetInfoFile(logger, NULL), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_INFO, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_INFO), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  ASSERT_EQ(logger, nullptr);

  ASSERT_EQ(std::fclose(fp), 0);
  EXPECT_EQ(CountLines(ReadFile(infofile)), 1);
  (void)std::remove(infofile.c_str());
#endif
}

TEST(SUNLoggerTest, SetDebugFileAndDisableWithNull)
{
#if SUNDIALS_LOGGING_LEVEL < SUNDIALS_LOGGING_DEBUG
  GTEST_SKIP() << "Debug logging not enabled in this build";
#else
  const std::string debugfile = "test_sundials_logger.debug_fp.debug";
  (void)std::remove(debugfile.c_str());

  FILE* fp = std::fopen(debugfile.c_str(), "w+");
  ASSERT_NE(fp, nullptr);

  SUNLogger logger = NULL;
  ASSERT_EQ(SUNLogger_Create(SUN_COMM_NULL, 0, &logger), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetDebugFile(logger, fp), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, "scope", "label",
                               "first"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_DEBUG), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_SetDebugFile(logger, NULL), SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, "scope", "label",
                               "second"),
            SUN_SUCCESS);
  ASSERT_EQ(SUNLogger_Flush(logger, SUN_LOGLEVEL_DEBUG), SUN_SUCCESS);

  ASSERT_EQ(SUNLogger_Destroy(&logger), SUN_SUCCESS);
  ASSERT_EQ(logger, nullptr);

  ASSERT_EQ(std::fclose(fp), 0);
  EXPECT_EQ(CountLines(ReadFile(debugfile)), 1);
  (void)std::remove(debugfile.c_str());
#endif
}
