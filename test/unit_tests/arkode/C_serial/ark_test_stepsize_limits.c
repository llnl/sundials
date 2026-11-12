/* -----------------------------------------------------------------------------
 * Regression test for applying updated absolute step-size limits at trial
 * selection.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(SUN_RCONST(0.0), ydot);
  return 0;
}

static int check(int retval, const char* name)
{
  if (retval == ARK_SUCCESS) { return 0; }
  fprintf(stderr, "%s returned %d\n", name, retval);
  return 1;
}

static int check_step(void* mem, N_Vector y, sunrealtype tout,
                      sunrealtype expected_h, sunrealtype expected_t)
{
  sunrealtype tret, hlast;
  if (check(ARKodeEvolve(mem, tout, y, &tret, ARK_ONE_STEP), "ARKodeEvolve") ||
      check(ARKodeGetLastStep(mem, &hlast), "ARKodeGetLastStep"))
  {
    return 1;
  }

  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF *
                          SUNMAX(SUN_RCONST(1.0), SUNRabs(expected_h));
  if (SUNRabs(hlast - expected_h) > tol ||
      SUNRabs((tret - expected_t) - expected_h) > tol)
  {
    fprintf(stderr,
            "Expected step " SUN_FORMAT_G " at t=" SUN_FORMAT_G
            ", got " SUN_FORMAT_G " at t=" SUN_FORMAT_G "\n",
            expected_h, expected_t, hlast, tret);
    return 1;
  }
  return 0;
}

static int check_underflowed_proposal(sunrealtype direction, SUNContext ctx,
                                      N_Vector y)
{
  void* mem = ERKStepCreate(rhs, SUN_RCONST(0.0), y, ctx);
  if (mem == NULL) { return 1; }

  int fail                    = 0;
  const sunrealtype initial_h = direction * SUN_RCONST(1.0);
  fail |= check(ARKodeSetInitStep(mem, initial_h), "ARKodeSetInitStep");
  fail |= check(ARKodeSetMinStep(mem, SUN_RCONST(0.25)), "ARKodeSetMinStep");
  fail |= check(ARKodeSetMaxStep(mem, SUN_RCONST(0.0)), "ARKodeSetMaxStep");
  fail |= check_step(mem, y, direction * SUN_RCONST(10.0), initial_h,
                     SUN_RCONST(0.0));

  /* Model a consumed proposal whose multiplication underflowed to signed zero. */
  ARKodeMem ark_mem = (ARKodeMem)mem;
  ark_mem->h        = initial_h;
  ark_mem->hprime   = SUN_RCONST(2.0) * initial_h;
  ark_mem->eta      = SUN_RCONST(0.0);

  const sunrealtype expected_h = direction * SUN_RCONST(0.25);
  fail |= check_step(mem, y, direction * SUN_RCONST(10.0), expected_h, initial_h);

  ARKodeFree(&mem);
  return fail;
}

int main(void)
{
  SUNContext ctx;
  if (SUNContext_Create(SUN_COMM_NULL, &ctx) != 0) { return 1; }

  N_Vector y = N_VNew_Serial(1, ctx);
  if (y == NULL)
  {
    SUNContext_Free(&ctx);
    return 1;
  }
  N_VConst(SUN_RCONST(1.0), y);

  void* mem = ERKStepCreate(rhs, SUN_RCONST(0.0), y, ctx);
  if (mem == NULL)
  {
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    return 1;
  }

  int fail = 0;
  fail |= check(ARKodeSetInitStep(mem, SUN_RCONST(1.0)), "ARKodeSetInitStep");
  fail |= check(ARKodeSetMaxFirstGrowth(mem, SUN_RCONST(2.0)),
                "ARKodeSetMaxFirstGrowth");
  fail |= check(ARKodeSetMaxGrowth(mem, SUN_RCONST(2.0)), "ARKodeSetMaxGrowth");
  fail |= check(ARKodeSetMinStep(mem, SUN_RCONST(0.0)), "ARKodeSetMinStep");
  fail |= check(ARKodeSetMaxStep(mem, SUN_RCONST(0.0)), "ARKodeSetMaxStep");
  fail |= check(ARKodeSetFixedStepBounds(mem, SUN_RCONST(0.0), SUN_RCONST(3.0)),
                "ARKodeSetFixedStepBounds");

  /* The first accepted step establishes an unchanged pending proposal. */
  fail |= check_step(mem, y, SUN_RCONST(10.0), SUN_RCONST(1.0), SUN_RCONST(0.0));
  ARKodeMem ark_mem = (ARKodeMem)mem;
  if (ark_mem->hprime != ark_mem->h)
  {
    fprintf(stderr,
            "Expected unchanged pending proposal before bound update\n");
    fail = 1;
  }

  /* Lower the maximum below the pending proposal: next trial must be 0.5. */
  fail |= check(ARKodeSetMaxStep(mem, SUN_RCONST(0.5)), "ARKodeSetMaxStep");
  fail |= check_step(mem, y, SUN_RCONST(10.0), SUN_RCONST(0.5), SUN_RCONST(1.0));

  /* Raise the minimum above the next pending proposal: next trial must be 1.5. */
  fail |= check(ARKodeSetMaxStep(mem, SUN_RCONST(2.0)), "ARKodeSetMaxStep");
  fail |= check(ARKodeSetMinStep(mem, SUN_RCONST(1.5)), "ARKodeSetMinStep");
  fail |= check_step(mem, y, SUN_RCONST(10.0), SUN_RCONST(1.5), SUN_RCONST(1.5));

  fail |= check_underflowed_proposal(SUN_RCONST(1.0), ctx, y);
  fail |= check_underflowed_proposal(SUN_RCONST(-1.0), ctx, y);

  ARKodeFree(&mem);
  N_VDestroy(y);
  SUNContext_Free(&ctx);
  return fail;
}
