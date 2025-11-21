/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * This is a 'private' header file and should not be used in user
 * code. It is subject to change without warning.
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * -----------------------------------------------------------------
 * Contains all error checking macros and private error handling API.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_ERRORS_IMPL_H
#define _SUNDIALS_ERRORS_IMPL_H

#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_export.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------------------------------------
 * Macros used in error handling
 * ---------------------------------------------------------------------------*/

/* ------------------------------------------------------------------
 * SUNDIALS __builtin_expect related macros.
 * These macros provide hints to the compiler that the condition
 * is typically false (or true) which may allow the compiler to
 * optimize.
 * -----------------------------------------------------------------*/

/* Hint to the compiler that the branch is unlikely to be taken */
#ifdef SUNDIALS_C_COMPILER_HAS_BUILTIN_EXPECT
#define SUNHintFalse(cond) __builtin_expect((cond), 0)
#else
#define SUNHintFalse(cond) (cond)
#endif

/* Hint to the compiler that the branch is likely to be taken */
#ifdef SUNDIALS_C_COMPILER_HAS_BUILTIN_EXPECT
#define SUNHintTrue(cond) __builtin_expect((cond), 1)
#else
#define SUNHintTrue(cond) (cond)
#endif

/* ------------------------------------------------------------------
 * SUNAssume
 *
 * This macro tells the compiler that the condition should be assumed
 * to be true. The consequence is that what happens if the assumption
 * is violated is undefined. If there is not compiler support for
 * assumptions, then we dont do anything as there is no reliable
 * way to avoid the condition being executed in all cases (such as
 * the condition being an opaque function call, which we have a lot of).
 * -----------------------------------------------------------------*/

#if defined(__cplusplus) && __cplusplus >= 202302L
#define SUNAssume(...) [[assume(__VA_ARGS__)]]
#elif defined(SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_ASSUME)
#define SUNAssume(...) __attribute__((assume(__VA_ARGS__)))
#elif defined(SUNDIALS_C_COMPILER_HAS_BUILTIN_ASSUME)
#define SUNAssume(...) __builtin_assume(__VA_ARGS__)
#elif defined(SUNDIALS_C_COMPILER_HAS_ASSUME)
#define SUNAssume(...) __assume(__VA_ARGS__)
#else
#define SUNAssume(...)
#endif

/* ----------------------------------------------------------------------------
 * SUNErrHandler_ definition.
 * ---------------------------------------------------------------------------*/

struct SUNErrHandler_
{
  SUNErrHandler previous; /* next error handler to call (singly linked-list) */
  SUNErrHandlerFn call;
  void* data;
};

/*
  This function creates a new SUNErrHandler object which is a node in a
  singly linked-lis of error handlers.

  :param eh_fn: An error handler callback function
  :param eh_data: A pointer that will be passed back to the error handler
  callback function
  :param eh_out: The new SUNErrHandler object

  :return: A SUNErrCode indicating success or failure
*/
SUNDIALS_EXPORT
SUNErrCode SUNErrHandler_Create(SUNErrHandlerFn eh_fn, void* eh_data,
                                SUNErrHandler* eh_out);

/*
  This function destroys and frees the memory for the given SUNErrHandler
  object.

  :param eh: A pointer to a SUNErrHandler object

  :return: void
*/
SUNDIALS_EXPORT
void SUNErrHandler_Destroy(SUNErrHandler* eh);

/*
  This function will print an error out to stderr. It is used as a fallback
  when the SUNContext is NULL or corrupt. It should not be used otherwise.

  :param line: the line number of the error
  :param func: the function in which the error occurred
  :param file: the file in which the error occurred
  :param msg: a message associated with the error
  :param args: the arguments to be provided to the format message

  :return: void
*/
SUNDIALS_EXPORT
void SUNGlobalFallbackErrHandler(int line, const char* func, const char* file,
                                 const char* msgfmt, SUNErrCode code, ...);

/*
  This function calls the error handlers registered with the SUNContext
  with the provided message.

  :param line: the line number of the error
  :param func: the function in which the error occurred
  :param file: the file in which the error occurred
  :param msgfmt: a message associated with the error with formatting
  :param code: the SUNErrCode for the error
  :param sunctx: a valid SUNContext object

  :return: void
*/
static inline void SUNHandleErrWithMsg(int line, const char* func,
                                       const char* file, const char* msg,
                                       SUNErrCode code, SUNContext sunctx)
{
  if (!sunctx) { SUNGlobalFallbackErrHandler(line, func, file, msg, code); }

  sunctx->last_err = code;
  SUNErrHandler eh = sunctx->err_handler;
  while (eh != NULL)
  {
    eh->call(line, func, file, msg, code, eh->data, sunctx);
    eh = eh->previous;
  }
}

/*
  This function calls the error handlers registered with the SUNContext
  with the provided format message.

  :param line: the line number of the error
  :param func: the function in which the error occurred
  :param file: the file in which the error occurred
  :param msgfmt: a message associated with the error with formatting
  :param code: the SUNErrCode for the error
  :param sunctx: a valid SUNContext object
  :param args: the arguments to be provided to the format message

  :return: void
*/
static inline void SUNHandleErrWithFmtMsg(int line, const char* func,
                                          const char* file, const char* msgfmt,
                                          SUNErrCode code, SUNContext sunctx, ...)
{
  size_t msglen;
  char* msg;
  va_list values;
  va_start(values, sunctx);
  msglen = (size_t)vsnprintf(NULL, (size_t)0, msgfmt, values); /* determine size
                                                                  of buffer
                                                                  needed */
  va_end(values);
  msg = (char*)malloc(msglen + 1);
  va_start(values, sunctx);
  vsnprintf(msg, msglen + 1, msgfmt, values);
  va_end(values);
  SUNHandleErrWithMsg(line, func, file, msg, code, sunctx);
  free(msg);
}

/*
  The SUNCTX_ macro expands to the name of the local SUNContext object
  defined by SUNFunctionBegin. SUNCTX_ should be used to reference the
  SUNContext inside of SUNDIALS functions.
 */
#define SUNCTX_ sunctx_local_scope_

/*
  The SUNFunctionBegin macro is used to declare the local SUNContext object to
  be used a function. It should be used at the start of every SUNDIALS
  functions.

  :param sunctx: the SUNContext to set the local SUNContext variable to
 */
#define SUNFunctionBegin(sunctx) \
  SUNContext SUNCTX_ = sunctx;   \
  (void)SUNCTX_

/* ----------------------------------------------------------------------------
 * SUNCheck* family of error checking macros
 *
 * We define several different versions of SUNCheck* macros to cover various
 * programming scenarios.
 *
 * SUNCheckCall* macros are used to check SUNDIALS function calls that
 * return a SUNErrCode.
 *
 * SUNCheckLastErr* macros are used to check SUNDIALS function calls that
 * do not return a SUNErrCode.
 * ---------------------------------------------------------------------------*/

/*
  SUNCheck evaluates the given expression and calls the error handler if it is
  false. It should be used to check expressions that do not imply stringent
  assumptions. If the expression should be strictly assumed as true, then use
  SUNAssert instead.

  Use SUNAssert macros to check for conditions that do not make sense e.g.,
  to check if malloc returned NULL.

  :param expr: an expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheck(expr, code)                                              \
  do {                                                                    \
    if (SUNHintFalse(!(expr)))                                            \
    {                                                                     \
      SUNHandleErrWithFmtMsg(__LINE__, __func__, __FILE__, "expected %s", \
                             code, SUNCTX_, #expr);                       \
      return code;                                                        \
    }                                                                     \
  }                                                                       \
  while (0)
#else
#define SUNCheck(expr, code)
#endif

/*
  SUNCheckNoRet is the same as SUNCheck but *does not return from the caller*.
  Use SUNAssert macros to check for conditions that do not make sense e.g.,
  to check if malloc returned NULL.

  :param expr: an expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckNoRet(expr, code)                                         \
  do {                                                                    \
    if (SUNHintFalse(!(expr)))                                            \
    {                                                                     \
      SUNHandleErrWithFmtMsg(__LINE__, __func__, __FILE__, "expected %s", \
                             code, SUNCTX_, #expr);                       \
    }                                                                     \
  }                                                                       \
  while (0)
#else
#define SUNCheckNoRet(expr, code)
#endif

/*
  SUNCheckNull is the same as SUNCheck but *returns NULL from the caller*.
  Use SUNAssert macros to check for conditions that do not make sense e.g.,
  to check if malloc returned NULL.

  :param expr: an expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckNull(expr, code)                                          \
  do {                                                                    \
    if (SUNHintFalse(!(expr)))                                            \
    {                                                                     \
      SUNHandleErrWithFmtMsg(__LINE__, __func__, __FILE__, "expected %s", \
                             code, SUNCTX_, #expr);                       \
      return NULL;                                                        \
    }                                                                     \
  }                                                                       \
  while (0)
#else
#define SUNCheckNull(expr, code)
#endif

/*
  SUNCheckNull is the same as SUNCheck but *returns void from the caller*.
  Use SUNAssert macros to check for conditions that do not make sense e.g.,
  to check if malloc returned NULL.

  :param expr: an expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckVoid(expr, code)                                          \
  do {                                                                    \
    if (SUNHintFalse(!(expr)))                                            \
    {                                                                     \
      SUNHandleErrWithFmtMsg(__LINE__, __func__, __FILE__, "expected %s", \
                             code, SUNCTX_, #expr);                       \
      return;                                                             \
    }                                                                     \
  }                                                                       \
  while (0)
#else
#define SUNCheckVoid(expr, code)
#endif

/*
   SUNCheckCallMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, call the error handler, **and then return the error code**.

   :param call: the function call
   :param msg: an error message
*/
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallMsg(call, msg)                           \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      return sun_chk_call_err_code_;                         \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallMsg(call, msg) (void)call
#endif

/*
   SUNCheckCallNoRetMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, and call the error handler. **It does not return**.

   :param call: the function call
   :param msg: an error message
*/
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallNoRetMsg(call, msg)                      \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallNoRetMsg(call, msg) (void)call
#endif

/*
   SUNCheckCallNullMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, call the error handler, **and then returns NULL**.

   :param call: the function call
   :param msg: an error message
*/
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallNullMsg(call, msg)                       \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      return NULL;                                           \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallNullMsg(call, msg) (void)call
#endif

/*
   SUNCheckCallNullMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, call the error handler, **and then returns void**.

   :param call: the function call
   :param msg: an error message
*/
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckCallVoidMsg(call, msg)                       \
  do {                                                       \
    SUNErrCode sun_chk_call_err_code_ = call;                \
    if (SUNHintFalse(sun_chk_call_err_code_ < 0))            \
    {                                                        \
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__, msg, \
                          sun_chk_call_err_code_, SUNCTX_);  \
      return;                                                \
    }                                                        \
  }                                                          \
  while (0)
#else
#define SUNCheckCallVoidMsg(call, msg) (void)call
#endif

/* These versions of SUNCheckCall do not take a custom message so a
   default message associated with the error code will be used. */
#define SUNCheckCall(call)      SUNCheckCallMsg(call, NULL)
#define SUNCheckCallNoRet(call) SUNCheckCallNoRetMsg(call, NULL)
#define SUNCheckCallNull(call)  SUNCheckCallNullMsg(call, NULL)
#define SUNCheckCallVoid(call)  SUNCheckCallVoidMsg(call, NULL)

/* SUNCheckLastErrMsg checks the last_err value in the SUNContext.
   If an error occurred, then it will log the error, set the last_err
   value, and call the error handler, **and then returns the code**. */
#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNCheckLastErrMsg(msg)                              \
  do {                                                       \
    SUNCheckCallMsg(SUNContext_PeekLastError(SUNCTX_), msg); \
  }                                                          \
  while (0)

/*
   SUNCheckLastErrNoRetMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, call the error handler. **It does not return.**

   :param msg: an error message
*/
#define SUNCheckLastErrNoRetMsg(msg)                              \
  do {                                                            \
    SUNCheckCallNoRetMsg(SUNContext_PeekLastError(SUNCTX_), msg); \
  }                                                               \
  while (0)

/*
   SUNCheckLastErrNullMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, call the error handler, **and then returns NULL**.

   :param msg: an error message
*/
#define SUNCheckLastErrNullMsg(msg)                              \
  do {                                                           \
    SUNCheckCallNullMsg(SUNContext_PeekLastError(SUNCTX_), msg); \
  }                                                              \
  while (0)

/*
   SUNCheckLastErrVoidMsg performs the SUNDIALS function call, and checks the
   returned error code. If an error occurred, then it will log the error, set the
   last_err value, call the error handler, **and then returns void**.

   :param msg: an error message
*/
#define SUNCheckLastErrVoidMsg(msg)                              \
  do {                                                           \
    SUNCheckCallVoidMsg(SUNContext_PeekLastError(SUNCTX_), msg); \
  }                                                              \
  while (0)
#else
#define SUNCheckLastErrNoRetMsg(msg)
#define SUNCheckLastErrMsg(msg)
#define SUNCheckLastErrVoidMsg(msg)
#define SUNCheckLastErrNullMsg(msg)
#endif

/* These versions of SUNCheckLastErr do not take a custom message so the
   default message associated with the error code will be used. */
#define SUNCheckLastErr()      SUNCheckLastErrMsg(NULL)
#define SUNCheckLastErrNoRet() SUNCheckLastErrNoRetMsg(NULL)
#define SUNCheckLastErrVoid()  SUNCheckLastErrVoidMsg(NULL)
#define SUNCheckLastErrNull()  SUNCheckLastErrNullMsg(NULL)

/*
  SUNAssert checks if an expression is true. It expands to SUNCheck in debug
  builds othewrise we try to expand it to an assumption, if the compiler
  supports assumptions, so that the compiler can make optimizations based on the
  assumption.

  :param expr: a expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNAssert(expr, code) SUNCheck(expr, code)
#else
#define SUNAssert(expr, code)
#endif

/*
  SUNAssertNoRet is the same as SUNAssert but it does not return from the
  caller.

  :param expr: a expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNAssertNoRet(expr, code) SUNCheckNoRet(expr, code)
#else
#define SUNAssertNoRet(expr, code)
#endif

/*
  SUNAssertNull is the same as SUNAssert but it *returns NULL from the caller*.

  :param expr: a expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNAssertNull(expr, code) SUNCheckNull(expr, code)
#else
#define SUNAssertNull(expr, code)
#endif

/*
  SUNAssertVoid is the same as SUNAssert but it *returns void from the caller*.

  :param expr: a expression to evaluate as true or false
  :param code: the error code to pass to the error handler if the expression is
  false
*/

#if defined(SUNDIALS_ENABLE_ERROR_CHECKS)
#define SUNAssertVoid(expr, code) SUNCheckVoid(expr, code)
#else
#define SUNAssertVoid(expr, code)
#endif

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_ERRORS_IMPL_H */
