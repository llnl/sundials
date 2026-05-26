#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>
#include <arkode/arkode_pdaestep.h>
#include <nvector/nvector_serial.h>

#define PARTITIONS 2

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SIN(x)  (sin(x))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x)  (sinf(x))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x)  (sinl(x))
#endif

static int f(sunrealtype t, N_Vector y, N_Vector w, N_Vector yp, N_Vector res, void* user_data)
{
  return 0;
}

static int h(sunrealtype t, N_Vector y, N_Vector w, N_Vector res, void* user_data)
{
  return 0;
}

int main(void) {
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  N_Vector x[] = {N_VNew_Serial(1, sunctx), N_VNew_Serial(1, sunctx)};
  N_VConst(SUN_RCONST(2.0), x[0]);
  N_VConst(SUN_RCONST(1.0), x[1]);
  N_Vector z[] = {N_VNew_Serial(1, sunctx), N_VNew_Serial(1, sunctx)};
  N_VConst(SUN_RCONST(1.0), z[0]);
  N_VConst(SUN_RCONST(2.0), z[1]);
  N_Vector w = N_VNew_Serial(1, sunctx);
  N_VConst(0.0, w);

  const sunrealtype s1 = SIN(SUN_RCONST(1.0));
  const sunrealtype s2 = SIN(SUN_RCONST(2.0));
  N_Vector xp[] = {N_VNew_Serial(1, sunctx), N_VNew_Serial(1, sunctx)};
  N_VConst(-SUN_RCONST(2.0) - s1, xp[0]);
  N_VConst(SUN_RCONST(2.0) * (SUN_RCONST(2.0) + s1), xp[1]);
  N_Vector zp[] = {N_VNew_Serial(1, sunctx), N_VNew_Serial(1, sunctx)};
  N_VConst(-SUN_RCONST(1.0) - s2, zp[0]);
  N_VConst((SUN_RCONST(1.0) + s2) / SUN_RCONST(2.0), zp[1]);
  N_Vector wp = N_VNew_Serial(1, sunctx);
  N_VConst(SUN_RCONST(1.5) * (SUN_RCONST(5.0) + SUN_RCONST(2.0) * s1 + s2), wp);

  N_Vector y = PDAEStepManyVector(x, z, w, PARTITIONS);
  N_Vector yp = PDAEStepManyVector(xp, zp, wp, PARTITIONS);

  PDAEStepComponentResFn component_res_fns[] = {f, f};
  PDAEStepAlgebraicResFn algebraic_res_fn = h;
  void *arkode_mem = PDAEStepCreate(component_res_fns, algebraic_res_fn, 0.0, y, yp, PARTITIONS, sunctx);

  for (int i = 0; i < PARTITIONS; i++) {
    N_VDestroy(x[i]);
    N_VDestroy(z[i]);
    N_VDestroy(xp[i]);
    N_VDestroy(zp[i]);
  }
  N_VDestroy(w);
  N_VDestroy(wp);
  N_VDestroy(y);
  N_VDestroy(yp);

  ARKodeFree(&arkode_mem);

  SUNContext_Free(&sunctx);
}
