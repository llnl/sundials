#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>
#include <arkode/arkode_pdaestep.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_manyvector.h>

#define PARTITIONS 2

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SIN(x)  (sin(x))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x)  (sinf(x))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x)  (sinl(x))
#endif

static sunrealtype f(sunrealtype x, sunrealtype z, sunrealtype w, sunrealtype xp) {
  return SIN(w - z) - x - xp;
}

static sunrealtype g(sunrealtype x, sunrealtype z, sunrealtype w) {
  return x * x + z * z + w * w - SUN_RCONST(5.0);
}

static int component_res(sunrealtype t, N_Vector y_vec, N_Vector w_vec, N_Vector yp_vec,
  N_Vector res_vec, void* user_data)
{
  // TODO(SBR): make helper functions to extract components
  sunrealtype x = N_VGetArrayPointer(N_VGetSubvector_ManyVector(y_vec, 0))[0];
  sunrealtype z = N_VGetArrayPointer(N_VGetSubvector_ManyVector(y_vec, 1))[0];
  sunrealtype xp = N_VGetArrayPointer(N_VGetSubvector_ManyVector(yp_vec, 0))[0];
  sunrealtype w = N_VGetArrayPointer(w_vec)[0];

  sunrealtype *res = N_VGetArrayPointer(res_vec);

  res[0] = f(x, z, w, xp);
  res[1] = g(x, z, w);

  return 0;
}

static sunrealtype h(sunrealtype x1, sunrealtype x2, sunrealtype z1, sunrealtype z2, sunrealtype w)
{
  return x1 - z1 + x2 - z2 + w;
}

static int algebraic_res(sunrealtype t, N_Vector y_vec, N_Vector w_vec, N_Vector res_vec, void* user_data)
{
  sunrealtype x1 = N_VGetArrayPointer(PDAEStepGetDifferentialSubvector(y_vec, 0))[0];
  sunrealtype z1 = N_VGetArrayPointer(PDAEStepGetAlgebraicSubvector(y_vec, 0))[0];
  sunrealtype x2 = N_VGetArrayPointer(PDAEStepGetDifferentialSubvector(y_vec, 1))[0];
  sunrealtype z2 = N_VGetArrayPointer(PDAEStepGetAlgebraicSubvector(y_vec, 1))[0];
  sunrealtype w = N_VGetArrayPointer(w_vec)[0];

  sunrealtype *res = N_VGetArrayPointer(res_vec);
  res[0] = g(x1, z1, w);
  res[1] = g(x1, z1, w);
  res[2] = h(x1, x2, z1, z2, w);

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
  const sunrealtype tspan[] = {SUN_RCONST(0.0), SUN_RCONST(1.0)};

  PDAEStepComponentResFn component_res_fns[] = {component_res, component_res};
  PDAEStepAlgebraicResFn algebraic_res_fn = algebraic_res;
  void *arkode_mem = PDAEStepCreate(component_res_fns, algebraic_res_fn, tspan[0], y, yp, PARTITIONS, sunctx);
  ARKodeSetFixedStep(arkode_mem, SUN_RCONST(0.01));

  sunrealtype tret;
  ARKodeEvolve(arkode_mem, tspan[1], y, &tret, ARK_NORMAL);

  N_VPrint(y);

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
