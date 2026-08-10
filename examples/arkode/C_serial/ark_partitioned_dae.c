#include <arkode/arkode.h>
#include <arkode/arkode_pdaestep.h>
#include <ida/ida.h>
#include <math.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#define PARTITIONS 2

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SIN(x) (sin(x))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x) (sinf(x))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x) (sinl(x))
#endif

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Private function to compute the max-norm solution error */
static void solution_error(N_Vector y, sunrealtype* err, sunrealtype* diff_err,
                           sunrealtype* alg_err);

static sunrealtype f(sunrealtype x, sunrealtype z, sunrealtype w, sunrealtype xp)
{
  return SIN(w - z) - x - xp;
}

static sunrealtype g(sunrealtype x, sunrealtype z, sunrealtype w)
{
  return x * x + z * z + w * w - SUN_RCONST(5.0);
}

static int component_res(int partition, sunrealtype t, N_Vector y_vec,
                         N_Vector w_vec, N_Vector yp_vec, N_Vector res_vec,
                         void* user_data)
{
  (void)partition;
  // TODO(SBR): make helper functions to extract components
  sunrealtype x  = N_VGetArrayPointer(N_VGetSubvector_ManyVector(y_vec, 0))[0];
  sunrealtype z  = N_VGetArrayPointer(N_VGetSubvector_ManyVector(y_vec, 1))[0];
  sunrealtype xp = N_VGetArrayPointer(N_VGetSubvector_ManyVector(yp_vec, 0))[0];
  sunrealtype w  = N_VGetArrayPointer(w_vec)[0];

  sunrealtype* x_res = N_VGetArrayPointer(N_VGetSubvector_ManyVector(res_vec, 0));
  sunrealtype* z_res = N_VGetArrayPointer(N_VGetSubvector_ManyVector(res_vec, 1));

  *x_res = f(x, z, w, xp);
  *z_res = g(x, z, w);

  return 0;
}

static sunrealtype h(sunrealtype x1, sunrealtype x2, sunrealtype z1,
                     sunrealtype z2, sunrealtype w)
{
  return x1 - z1 + x2 - z2 + w;
}

static int algebraic_res(sunrealtype t, N_Vector y_vec, N_Vector w_vec,
                         N_Vector res_vec, void* user_data)
{
  sunrealtype x1 =
    N_VGetArrayPointer(PDAEStepGetDifferentialSubvector(y_vec, 0))[0];
  sunrealtype z1 = N_VGetArrayPointer(PDAEStepGetAlgebraicSubvector(y_vec, 0))[0];
  sunrealtype x2 =
    N_VGetArrayPointer(PDAEStepGetDifferentialSubvector(y_vec, 1))[0];
  sunrealtype z2 = N_VGetArrayPointer(PDAEStepGetAlgebraicSubvector(y_vec, 1))[0];
  sunrealtype w = N_VGetArrayPointer(w_vec)[0];

  // TODO(SBR): make helper functions to extract components
  sunrealtype* z1_res = N_VGetArrayPointer(N_VGetSubvector_ManyVector(res_vec, 0));
  sunrealtype* z2_res = N_VGetArrayPointer(N_VGetSubvector_ManyVector(res_vec, 1));
  sunrealtype* w_res = N_VGetArrayPointer(N_VGetSubvector_ManyVector(res_vec, 2));
  *z1_res = g(x1, z1, w);
  *z2_res = g(x2, z2, w);
  *w_res  = h(x1, x2, z1, z2, w);

  return 0;
}

int main(void)
{
  int flag                                   = 0;
  SUNContext sunctx                          = NULL;
  N_Vector x[PARTITIONS]                     = {NULL};
  N_Vector z[PARTITIONS]                     = {NULL};
  N_Vector xp[PARTITIONS]                    = {NULL};
  N_Vector zp[PARTITIONS]                    = {NULL};
  N_Vector w                                 = NULL;
  N_Vector wp                                = NULL;
  N_Vector y                                 = NULL;
  N_Vector yp                                = NULL;
  N_Vector alg_template                      = NULL;
  N_Vector alg_res                           = NULL;
  void* arkode_mem                           = NULL;
  SUNNonlinearSolver alg_nls                 = NULL;
  SUNLinearSolver alg_linear_solver          = NULL;
  SUNLinearSolver linear_solvers[PARTITIONS] = {NULL};
  sunrealtype tret                           = SUN_RCONST(0.0);

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&flag, "SUNContext_Create", 1)) { return 1; }

  x[0] = N_VNew_Serial(1, sunctx);
  x[1] = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)x[0], "N_VNew_Serial", 0)) { return 1; }
  if (check_retval((void*)x[1], "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(2.0), x[0]);
  N_VConst(SUN_RCONST(1.0), x[1]);
  z[0] = N_VNew_Serial(1, sunctx);
  z[1] = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)z[0], "N_VNew_Serial", 0)) { return 1; }
  if (check_retval((void*)z[1], "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(1.0), z[0]);
  N_VConst(SUN_RCONST(2.0), z[1]);
  w = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)w, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(0.0, w);

  const sunrealtype s1 = SIN(SUN_RCONST(1.0));
  const sunrealtype s2 = SIN(SUN_RCONST(2.0));
  xp[0]                = N_VNew_Serial(1, sunctx);
  xp[1]                = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)xp[0], "N_VNew_Serial", 0)) { return 1; }
  if (check_retval((void*)xp[1], "N_VNew_Serial", 0)) { return 1; }
  N_VConst(-SUN_RCONST(2.0) - s1, xp[0]);
  N_VConst(SUN_RCONST(2.0) * (SUN_RCONST(2.0) + s1), xp[1]);
  zp[0] = N_VNew_Serial(1, sunctx);
  zp[1] = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)zp[0], "N_VNew_Serial", 0)) { return 1; }
  if (check_retval((void*)zp[1], "N_VNew_Serial", 0)) { return 1; }
  N_VConst(-SUN_RCONST(1.0) - s2, zp[0]);
  N_VConst((SUN_RCONST(1.0) + s2) / SUN_RCONST(2.0), zp[1]);
  wp = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)wp, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(1.5) * (SUN_RCONST(5.0) + SUN_RCONST(2.0) * s1 + s2), wp);

  y = PDAEStepManyVector(x, z, w, PARTITIONS);
  if (check_retval((void*)y, "PDAEStepManyVector", 0)) { return 1; }
  yp = PDAEStepManyVector(xp, zp, wp, PARTITIONS);
  if (check_retval((void*)yp, "PDAEStepManyVector", 0)) { return 1; }
  const sunrealtype tspan[] = {SUN_RCONST(0.0), SUN_RCONST(1.0)};

  PDAEStepAlgebraicResFn algebraic_res_fn = algebraic_res;
  arkode_mem = PDAEStepCreate(component_res, algebraic_res_fn, tspan[0], y, yp,
                              PARTITIONS, sunctx);
  if (check_retval((void*)arkode_mem, "PDAEStepCreate", 0)) { return 1; }
  flag = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6));
  if (check_retval(&flag, "ARKodeSStolerances", 1)) { return 1; }
  flag = ARKodeSetFixedStep(arkode_mem, SUN_RCONST(0.01));
  if (check_retval(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = PDAEStepSetMaxNonlinIters(arkode_mem, 12);
  if (check_retval(&flag, "PDAEStepSetMaxNonlinIters", 1)) { return 1; }

  flag = PDAEStepGetAlgebraicVectorTemplate(arkode_mem, &alg_template);
  if (check_retval(&flag, "PDAEStepGetAlgebraicVectorTemplate", 1))
  {
    return 1;
  }

  alg_nls = SUNNonlinSol_Newton(alg_template, sunctx);
  if (check_retval((void*)alg_nls, "SUNNonlinSol_Newton", 0)) { return 1; }
  flag = PDAEStepSetNonlinearSolver(arkode_mem, alg_nls);
  if (check_retval(&flag, "PDAEStepSetNonlinearSolver", 1)) { return 1; }

  alg_linear_solver = SUNLinSol_SPGMR(alg_template, SUN_PREC_NONE, 0, sunctx);
  if (check_retval((void*)alg_linear_solver, "SUNLinSol_SPGMR", 0))
  {
    return 1;
  }
  flag = PDAEStepSetLinearSolver(arkode_mem, alg_linear_solver, NULL);
  if (check_retval(&flag, "PDAEStepSetLinearSolver", 1)) { return 1; }

  for (int i = 0; i < PARTITIONS; i++)
  {
    N_Vector linear_solver_temp = NULL;
    flag = PDAEStepGetPartitionVectorTemplate(arkode_mem, i, &linear_solver_temp);
    if (check_retval(&flag, "PDAEStepGetPartitionVectorTemplate", 1))
    {
      return 1;
    }
    linear_solvers[i] = SUNLinSol_SPGMR(linear_solver_temp, SUN_PREC_NONE, 0,
                                        sunctx);
    if (check_retval((void*)linear_solvers[i], "SUNLinSol_SPGMR", 0))
    {
      return 1;
    }
    void* ida_mem = NULL;
    flag          = PDAEStepGetPartitionIntegrator(arkode_mem, i, &ida_mem);
    if (check_retval(&flag, "PDAEStepGetPartitionIntegrator", 1)) { return 1; }
    flag = IDASetLinearSolver(ida_mem, linear_solvers[i], NULL);
    if (check_retval(&flag, "IDASetLinearSolver", 1)) { return 1; }
  }

  flag = ARKodeEvolve(arkode_mem, tspan[1], y, &tret, ARK_NORMAL);
  if (check_retval(&flag, "ARKodeEvolve", 1)) { return 1; }

  alg_res = N_VClone(alg_template);
  if (check_retval((void*)alg_res, "N_VClone", 0)) { return 1; }
  flag = algebraic_res(tret, y, w, alg_res, NULL);
  if (check_retval(&flag, "algebraic_res", 1)) { return 1; }
  sunrealtype alg_res_norm = N_VMaxNorm(alg_res);
  if (alg_res_norm > SUN_RCONST(1.0e-8))
  {
    fprintf(stderr, "Final algebraic residual norm is too large: %.16g\n",
            (double)alg_res_norm);
    flag = 1;
  }

  N_VPrint(y);

  sunrealtype err      = SUN_RCONST(0.0);
  sunrealtype diff_err = SUN_RCONST(0.0);
  sunrealtype alg_err  = SUN_RCONST(0.0);
  solution_error(y, &err, &diff_err, &alg_err);
  printf("\nSolution error at t = %.16g\n", (double)tret);
  printf("  Max-norm error      = %.16g\n", (double)err);
  printf("  Differential error  = %.16g\n", (double)diff_err);
  printf("  Algebraic error     = %.16g\n", (double)alg_err);

  ARKodeFree(&arkode_mem);
  if (alg_linear_solver) { SUNLinSolFree(alg_linear_solver); }
  if (alg_nls) { SUNNonlinSolFree(alg_nls); }
  if (alg_res) { N_VDestroy(alg_res); }
  for (int i = 0; i < PARTITIONS; i++)
  {
    N_VDestroy(x[i]);
    N_VDestroy(z[i]);
    N_VDestroy(xp[i]);
    N_VDestroy(zp[i]);
    if (linear_solvers[i]) { SUNLinSolFree(linear_solvers[i]); }
  }
  N_VDestroy(w);
  N_VDestroy(wp);
  N_VDestroy(y);
  N_VDestroy(yp);

  SUNContext_Free(&sunctx);
  return flag;
}

static void solution_error(N_Vector y, sunrealtype* err, sunrealtype* diff_err,
                           sunrealtype* alg_err)
{
  sunrealtype x0 = N_VGetArrayPointer(PDAEStepGetDifferentialSubvector(y, 0))[0];
  sunrealtype x1 = N_VGetArrayPointer(PDAEStepGetDifferentialSubvector(y, 1))[0];
  sunrealtype z0 = N_VGetArrayPointer(PDAEStepGetAlgebraicSubvector(y, 0))[0];
  sunrealtype z1 = N_VGetArrayPointer(PDAEStepGetAlgebraicSubvector(y, 1))[0];
  sunrealtype w  = N_VGetArrayPointer(PDAEStepGetCouplingSubvector(y))[0];

  const sunrealtype x0_exact = SUN_RCONST(0.81813251372736039768);
  const sunrealtype x1_exact = SUN_RCONST(0.23501133358643333385);
  const sunrealtype z0_exact = SUN_RCONST(1.2467607076859356555);
  const sunrealtype z1_exact = SUN_RCONST(1.4725904879949859286);
  const sunrealtype w_exact  = SUN_RCONST(1.6662073483671278526);

  *diff_err = SUNMAX(SUNRabs(x0 - x0_exact), SUNRabs(x1 - x1_exact));
  *alg_err  = SUNMAX(SUNRabs(w - w_exact), SUNRabs(z0 - z0_exact));
  *alg_err  = SUNMAX(*alg_err, SUNRabs(z1 - z1_exact));
  *err      = SUNMAX(*diff_err, *alg_err);
}

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
