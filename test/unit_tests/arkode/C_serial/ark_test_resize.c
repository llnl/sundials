/*-----------------------------------------------------------------
 * Regression tests for ARKodeResize initialization state.
 *---------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_impl.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_macros.h>

static ARKTimestepInitFn delegated_step_init;
static int step_init_calls;
static int last_init_type;

static int rhs(SUNDIALS_MAYBE_UNUSED sunrealtype t, N_Vector y, N_Vector ydot,
               SUNDIALS_MAYBE_UNUSED void* user_data)
{
  sunrealtype* y_data = N_VGetArrayPointer_Serial(y);
  sunrealtype* ydot_data = N_VGetArrayPointer_Serial(ydot);
  for (sunindextype i = 0; i < N_VGetLength(y); i++) { ydot_data[i] = -y_data[i]; }
  return 0;
}

static int step_init_spy(ARKodeMem ark_mem, int init_type)
{
  step_init_calls++;
  last_init_type = init_type;
  return delegated_step_init(ark_mem, init_type);
}

static int failing_resize(SUNDIALS_MAYBE_UNUSED ARKodeMem ark_mem,
                          SUNDIALS_MAYBE_UNUSED N_Vector ynew,
                          SUNDIALS_MAYBE_UNUSED sunrealtype hscale,
                          SUNDIALS_MAYBE_UNUSED sunrealtype t0,
                          SUNDIALS_MAYBE_UNUSED ARKVecResizeFn resize,
                          SUNDIALS_MAYBE_UNUSED void* resize_data)
{
  return ARK_MEM_FAIL;
}

static int check(int condition, const char* message)
{
  if (condition) { return 0; }
  fprintf(stderr, "Resize regression failed: %s\n", message);
  return 1;
}

static void install_step_init_spy(void* mem)
{
  ARKodeMem ark_mem   = (ARKodeMem)mem;
  delegated_step_init = ark_mem->step_init;
  ark_mem->step_init  = step_init_spy;
  step_init_calls     = 0;
  last_init_type      = -1;
}

static N_Vector make_vector(sunindextype length, SUNContext ctx)
{
  N_Vector y = N_VNew_Serial(length, ctx);
  if (y != NULL) { N_VConst(SUN_RCONST(1.0), y); }
  return y;
}

static int evolve(void* mem, N_Vector y, sunrealtype tout)
{
  sunrealtype tret;
  return ARKodeEvolve(mem, tout, y, &tret, ARK_NORMAL);
}

static int check_solution(N_Vector y, sunrealtype expected, const char* message)
{
  sunindextype i;
  sunrealtype tol = SUN_RCONST(1.0e-3);
  sunrealtype* y_data = N_VGetArrayPointer(y);

  for (i = 0; i < N_VGetLength(y); i++)
  {
    sunrealtype value = y_data[i];
    if (!isfinite(value) || SUNRabs(value - expected) > tol)
    {
      fprintf(stderr, "Resize regression failed: %s (component %ld)\n", message,
              (long)i);
      return 1;
    }
  }
  return 0;
}

static int preinitialization_resize_body(void* mem, N_Vector y1, N_Vector y2)
{
  int fail = 0;
  install_step_init_spy(mem);

  fail |= check(ARKodeReset(mem, SUN_RCONST(0.0), y1) == ARK_SUCCESS,
                "pre-initialization Reset");
  fail |= check(ARKodeResize(mem, y2, SUN_RCONST(1.0), SUN_RCONST(0.0), NULL,
                             NULL) == ARK_SUCCESS,
                "pre-initialization Resize");
  ARKodeMem ark_mem = (ARKodeMem)mem;
  fail |= check(!ark_mem->initialized && ark_mem->init_type == FIRST_INIT &&
                  !ark_mem->preallocated && ark_mem->initsetup &&
                  ark_mem->firststage,
                "pre-initialization Resize preserves FIRST_INIT");
  fail |= check(evolve(mem, y2, SUN_RCONST(0.1)) == ARK_SUCCESS,
                "evolve after pre-initialization Resize");
  fail |= check_solution(y2, SUNRexp(SUN_RCONST(-0.1)),
                         "solution after pre-initialization Resize");
  fail |= check(step_init_calls == 1 && last_init_type == FIRST_INIT,
                "first setup after pre-initialization Resize");
  return fail;
}

static int test_preinitialization_resize(SUNContext ctx)
{
  int fail    = 0;
  void* mem   = NULL;
  N_Vector y1 = make_vector(1, ctx);
  N_Vector y2 = make_vector(2, ctx);
  if (y1 == NULL || y2 == NULL) { fail = 1; }
  else
  {
    mem = ERKStepCreate(rhs, SUN_RCONST(0.0), y1, ctx);
    if (mem == NULL) { fail = 1; }
    else
    {
      fail |= preinitialization_resize_body(mem, y1, y2);
    }
  }

  if (mem != NULL) { ARKodeFree(&mem); }
  if (y1 != NULL) { N_VDestroy(y1); }
  if (y2 != NULL) { N_VDestroy(y2); }
  return fail;
}

static int preallocated_preinitialization_resize_body(void* mem, N_Vector y2)
{
  int fail = 0;
  install_step_init_spy(mem);

  fail |= check(ARKodeInit(mem) == ARK_SUCCESS, "ARKodeInit before Resize");
  fail |= check(ARKodeResize(mem, y2, SUN_RCONST(1.0), SUN_RCONST(0.0), NULL,
                             NULL) == ARK_SUCCESS,
                "preallocated pre-initialization Resize");
  ARKodeMem ark_mem = (ARKodeMem)mem;
  fail |= check(!ark_mem->initialized && ark_mem->init_type == FIRST_INIT &&
                  ark_mem->preallocated,
                "pre-initialization Resize preserves preallocation");
  fail |= check(evolve(mem, y2, SUN_RCONST(0.1)) == ARK_SUCCESS,
                "evolve after preallocated pre-initialization Resize");
  fail |=
    check_solution(y2, SUNRexp(SUN_RCONST(-0.1)),
                   "solution after preallocated pre-initialization Resize");
  return fail;
}

static int test_preallocated_preinitialization_resize(SUNContext ctx)
{
  int fail    = 0;
  void* mem   = NULL;
  N_Vector y1 = make_vector(1, ctx);
  N_Vector y2 = make_vector(2, ctx);
  if (y1 == NULL || y2 == NULL) { fail = 1; }
  else
  {
    mem = ERKStepCreate(rhs, SUN_RCONST(0.0), y1, ctx);
    if (mem == NULL) { fail = 1; }
    else
    {
      fail |= preallocated_preinitialization_resize_body(mem, y2);
    }
  }

  if (mem != NULL) { ARKodeFree(&mem); }
  if (y1 != NULL) { N_VDestroy(y1); }
  if (y2 != NULL) { N_VDestroy(y2); }
  return fail;
}

static int postinitialization_resize_body(void* mem, N_Vector y1, N_Vector y2,
                                          int preallocate, int reset_first)
{
  int fail = 0;
  install_step_init_spy(mem);

  if (preallocate)
  {
    fail |= check(ARKodeInit(mem) == ARK_SUCCESS, "ARKodeInit");
  }
  fail |= check(evolve(mem, y1, SUN_RCONST(0.1)) == ARK_SUCCESS,
                "initial evolve");
  if (reset_first)
  {
    fail |= check(ARKodeReset(mem, SUN_RCONST(0.1), y1) == ARK_SUCCESS,
                  "post-initialization Reset");
  }
  fail |= check(ARKodeResize(mem, y2, SUN_RCONST(1.0), SUN_RCONST(0.1), NULL,
                             NULL) == ARK_SUCCESS,
                "post-initialization Resize");
  ARKodeMem ark_mem = (ARKodeMem)mem;
  fail |= check(ark_mem->initialized && ark_mem->init_type == RESIZE_INIT &&
                  !ark_mem->preallocated && ark_mem->initsetup &&
                  ark_mem->firststage,
                "post-initialization Resize selects RESIZE_INIT");

  /* A second Resize must not replace the pending resize initialization. */
  fail |= check(ARKodeResize(mem, y2, SUN_RCONST(1.0), SUN_RCONST(0.1), NULL,
                             NULL) == ARK_SUCCESS,
                "repeated post-initialization Resize");
  fail |= check(ark_mem->initialized && ark_mem->init_type == RESIZE_INIT &&
                  !ark_mem->preallocated,
                "repeated Resize preserves RESIZE_INIT");
  fail |= check(evolve(mem, y2, SUN_RCONST(0.2)) == ARK_SUCCESS,
                "evolve after post-initialization Resize");
  fail |= check_solution(y2, SUNRexp(SUN_RCONST(-0.1)),
                         "solution after post-initialization Resize");
  fail |= check(step_init_calls == 2 && last_init_type == RESIZE_INIT,
                "next setup invokes step_init with RESIZE_INIT");
  return fail;
}

static int test_postinitialization_resize(SUNContext ctx, int preallocate,
                                          int reset_first)
{
  int fail    = 0;
  void* mem   = NULL;
  N_Vector y1 = make_vector(1, ctx);
  N_Vector y2 = make_vector(2, ctx);
  if (y1 == NULL || y2 == NULL) { fail = 1; }
  else
  {
    mem = ERKStepCreate(rhs, SUN_RCONST(0.0), y1, ctx);
    if (mem == NULL) { fail = 1; }
    else
    {
      fail |= postinitialization_resize_body(mem, y1, y2, preallocate,
                                             reset_first);
    }
  }

  if (mem != NULL) { ARKodeFree(&mem); }
  if (y1 != NULL) { N_VDestroy(y1); }
  if (y2 != NULL) { N_VDestroy(y2); }
  return fail;
}

/* Objects shared between the implicit test body and its cleanup */
typedef struct
{
  void* mem;
  SUNMatrix A1, A2;
  SUNLinearSolver LS1, LS2;
} implicit_resize_objects;

static int implicit_resize_body(implicit_resize_objects* obj, N_Vector y1,
                                N_Vector y2, SUNContext ctx)
{
  int fail = 0;

  obj->mem = ARKStepCreate(NULL, rhs, SUN_RCONST(0.0), y1, ctx);
  if (obj->mem == NULL) { return 1; }
  void* mem = obj->mem;

  fail |= check(ARKodeSetLinearSolver(mem, obj->LS1, obj->A1) == ARK_SUCCESS,
                "set initial implicit linear solver");
  install_step_init_spy(mem);
  fail |= check(ARKodeInit(mem) == ARK_SUCCESS,
                "preallocate implicit stepper before Resize");
  fail |= check(evolve(mem, y1, SUN_RCONST(0.1)) == ARK_SUCCESS,
                "initial implicit evolve");
  fail |= check(ARKodeResize(mem, y2, SUN_RCONST(1.0), SUN_RCONST(0.1), NULL,
                             NULL) == ARK_SUCCESS,
                "implicit post-initialization Resize");
  ARKodeMem ark_mem = (ARKodeMem)mem;
  fail |= check(ark_mem->init_type == RESIZE_INIT && !ark_mem->preallocated,
                "implicit Resize invalidates preallocation");

  /* Caller-owned matrix and linear solver must match the grown vector. */
  obj->A2  = SUNDenseMatrix(2, 2, ctx);
  obj->LS2 = SUNLinSol_Dense(y2, obj->A2, ctx);
  if (obj->A2 == NULL || obj->LS2 == NULL) { return 1; }
  fail |= check(ARKodeSetLinearSolver(mem, obj->LS2, obj->A2) == ARK_SUCCESS,
                "replace implicit linear solver after Resize");
  SUNLinSolFree(obj->LS1);
  SUNMatDestroy(obj->A1);
  obj->LS1 = NULL;
  obj->A1  = NULL;
  fail |= check(evolve(mem, y2, SUN_RCONST(0.2)) == ARK_SUCCESS,
                "implicit evolve at grown dimension");
  fail |= check_solution(y2, SUNRexp(SUN_RCONST(-0.1)),
                         "implicit solution after Resize");
  fail |= check(last_init_type == RESIZE_INIT,
                "implicit stepper receives RESIZE_INIT");
  return fail;
}

static int test_implicit_resize(SUNContext ctx)
{
  int fail                    = 0;
  N_Vector y1                 = make_vector(1, ctx);
  N_Vector y2                 = make_vector(2, ctx);
  implicit_resize_objects obj = {NULL, NULL, NULL, NULL, NULL};
  if (y1 == NULL || y2 == NULL) { fail = 1; }
  else
  {
    obj.A1  = SUNDenseMatrix(1, 1, ctx);
    obj.LS1 = SUNLinSol_Dense(y1, obj.A1, ctx);
    if (obj.A1 == NULL || obj.LS1 == NULL) { fail = 1; }
    else
    {
      fail |= implicit_resize_body(&obj, y1, y2, ctx);
    }
  }

  if (obj.mem != NULL) { ARKodeFree(&obj.mem); }
  if (obj.LS1 != NULL) { SUNLinSolFree(obj.LS1); }
  if (obj.LS2 != NULL) { SUNLinSolFree(obj.LS2); }
  if (obj.A1 != NULL) { SUNMatDestroy(obj.A1); }
  if (obj.A2 != NULL) { SUNMatDestroy(obj.A2); }
  if (y1 != NULL) { N_VDestroy(y1); }
  if (y2 != NULL) { N_VDestroy(y2); }
  return fail;
}

static int resize_failure_body(void* mem, N_Vector y1, N_Vector y2)
{
  int fail = 0;
  if (evolve(mem, y1, SUN_RCONST(0.1)) != ARK_SUCCESS) { return 1; }
  ARKodeMem ark_mem    = (ARKodeMem)mem;
  ark_mem->step_resize = failing_resize;
  fail |= check(ARKodeResize(mem, y2, SUN_RCONST(1.0), SUN_RCONST(0.1), NULL,
                             NULL) == ARK_MEM_FAIL,
                "stepper resize failure return code");
  fail |= check(ark_mem->initialized && ark_mem->init_type == RESIZE_INIT &&
                  !ark_mem->preallocated && ark_mem->initsetup &&
                  ark_mem->firststage,
                "state after stepper resize failure");
  return fail;
}

static int test_resize_failure(SUNContext ctx)
{
  int fail    = 0;
  void* mem   = NULL;
  N_Vector y1 = make_vector(1, ctx);
  N_Vector y2 = make_vector(2, ctx);
  if (y1 == NULL || y2 == NULL) { fail = 1; }
  else
  {
    mem = ERKStepCreate(rhs, SUN_RCONST(0.0), y1, ctx);
    if (mem == NULL) { fail = 1; }
    else
    {
      fail |= resize_failure_body(mem, y1, y2);
    }
  }

  if (mem != NULL) { ARKodeFree(&mem); }
  if (y1 != NULL) { N_VDestroy(y1); }
  if (y2 != NULL) { N_VDestroy(y2); }
  return fail;
}

int main(void)
{
  SUNContext ctx;
  int fail = 0;
  if (SUNContext_Create(SUN_COMM_NULL, &ctx) != SUN_SUCCESS) { return 1; }

  fail |= test_preinitialization_resize(ctx);
  fail |= test_preallocated_preinitialization_resize(ctx);
  fail |= test_postinitialization_resize(ctx, 0, 0);
  fail |= test_postinitialization_resize(ctx, 1, 0);
  fail |= test_postinitialization_resize(ctx, 0, 1);
  fail |= test_postinitialization_resize(ctx, 1, 1);
  fail |= test_implicit_resize(ctx);
  fail |= test_resize_failure(ctx);

  SUNContext_Free(&ctx);
  return fail;
}
