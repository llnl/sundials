/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 * The following is a simple example problem based off of ark_brusselator.c.
 *
 * We simulate a scenario where a set of independent ODEs are batched together
 * to form a larger system. Each independent ODE system has 3 components,
 * Y = [u,v,w], satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions Y0 = [u0,v0,w0].
 * The problem is stiff. We have 3 different testing scenarios:
 *
 * Reactor 0:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change
 *    during the first 0.2 time units, followed by a slow and
 *    smooth evolution.
 *
 * Reactor 1:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients
 *    during the first 0.3 time units, and all then proceed very
 *    smoothly for the remainder of the simulation.

 * Reactor 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5
 *    within a few steps.  All values proceed smoothly until
 *    around t=6.5, when both u and v undergo a sharp transition,
 *    with u increaseing from around 0.5 to 5 and v decreasing
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * This program solves the problem with the BDF method, Newton iteration, a
 * user-supplied Jacobian routine, and since the grouping of the independent
 * systems results in a block diagonal linear system, with the MAGMA
 * SUNLinearSolver which supports batched iterative methods. 100 outputs are
 * printed at equal intervals, and run statistics are printed at the end.
 *
 * The program takes three optional arguments, the number of independent ODE
 * systems (i.e. number of batches), the linear solver type
 * (MAGMA batched LU, non-batched GMRES with the Jacobian computed by
 * difference quotients, or non-batched GMRES with analytical Jacobian), and
 * the test type (uniform_1, uniform_2, uniform_3, or random).
 *
 *    ./cv_bruss_batched_magma [num_batches] [solver_type] [test_type]
 *
 * Options:
 *    num_batches <int>
 *    solver_type:
 *       0 - MAGMA batched GMRES
 *       1 - SUNDIALS non-batched GMRES with difference quotients Jacobian
 *    test_type:
 *       0 - all batches are the same reactor type 0
 *       1 - all batches are the same reactor type 1
 *       2 - all batches are the same reactor type 2
 *       3 - random distribution of reactors
 *       4 - all reactor type 0 but with random distribution of ep (stiffness)
 *       5 - all reactor type 1 but with random distribution of ep (stiffness)
 *       6 - all reactor type 2 but with random distribution of ep (stiffness)
 *       7 - random distribution of reactors and ep (stiffness)
 * --------------------------------------------------------------------------*/

#include <cstdio>
#include <cvode/cvode.h> /* prototypes for CVODE fcts., consts.           */
#include <memory>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_magmadense.h> /* access to MAGMA dense SUNLinearSolver         */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to GMRES SUNLinearSolver               */
#include <sunmatrix/sunmatrix_magmadense.h> /* access to the MAGMA dense SUNMatrix           */
#include <vector>

/* Convenience macro to help support both HIP and CUDA */
#if defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HIP_OR_CUDA(a, b) a
#elif defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HIP_OR_CUDA(a, b) b
#else
#define HIP_OR_CUDA(a, b) ((void)0);
#endif

/* Include the appropriate vector and memory helper based on the GPU target */
#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#include <nvector/nvector_hip.h>
#include <sunmemory/sunmemory_hip.h>
#else
#error SUNDIALS_MAGMA_BACKENDS value is not recognized
#endif

/* Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

__global__ static void f_kernel(sunrealtype t, sunrealtype* y, sunrealtype* ydot,
                                sunrealtype* A, sunrealtype* B, sunrealtype* Ep,
                                int nbatches, int batchSize);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

__global__ static void j_kernel(sunrealtype* ydata, sunrealtype* Jdata,
                                sunrealtype* A, sunrealtype* B, sunrealtype* Ep,
                                int nbatches, int batchSize);

/* Private function to output results */
static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3);

/* Private function to print final statistics */
static void PrintFinalStats(void* cvode_mem, SUNLinearSolver LS);

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Helper class for Arrays. Built on top of the SUNMemoryHelper. */
template<typename T, typename I>
class Array
{
public:
  Array(I size, SUNMemoryHelper helper) : helper_(helper), mem_(nullptr)
  {
    SUNMemoryHelper_Alloc(helper, &mem_, size * sizeof(T), SUNMEMTYPE_UVM, NULL);
  }

  Array(SUNMemory mem, SUNMemoryHelper helper) : helper_(helper), mem_(mem) {}

  T& operator[](int index) { return static_cast<T*>(mem_->ptr)[index]; }

  T* get() { return static_cast<T*>(mem_->ptr); }

private:
  SUNMemory mem_;
  SUNMemoryHelper helper_;
};

using RealArray = Array<sunrealtype, int>;

/* User data structure. This will be available
   in SUNDIALS callback functions. */
struct UserData
{
  UserData(int nbatches, int batchSize, SUNMemoryHelper h)
    : nbatches(nbatches),
      batchSize(batchSize),
      u0{nbatches, h},
      v0{nbatches, h},
      w0{nbatches, h},
      a{nbatches, h},
      b{nbatches, h},
      ep{nbatches, h}
  {}

  int nbatches;         /* number of chemical networks */
  int batchSize;        /* size of each network        */
  RealArray u0, v0, w0; /* initial conditions */
  RealArray a, b;       /* chemical concentrations that are constant */
  RealArray ep;
};

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char* argv[])
{
  const sunrealtype T0 = SUN_RCONST(0.0);  /* initial time                   */
  const sunrealtype Tf = SUN_RCONST(10.0); /* final time                    */
  const sunrealtype dTout = SUN_RCONST(1.0); /* time between outputs          */
  const int Nt = (int)ceil(Tf / dTout);      /* number of output times        */
  const sunrealtype reltol = 1.0e-6;         /* relative integrator tolerance */
  int retval;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;

  y = abstol = NULL;
  A          = NULL;
  LS         = NULL;
  cvode_mem  = NULL;

  /* Create the SUNDIALS context */
  sundials::Context sunctx;

  /* Create a SUNMemoryHelper for HIP or CUDA depending on the target.
     This will be used underneath the Arrays and the Magma solver. */
  SUNMemoryHelper memhelper =
    HIP_OR_CUDA(SUNMemoryHelper_Hip(sunctx);, SUNMemoryHelper_Cuda(sunctx);)
    /* Set defaults */
    int batchSize = 3;
  int nbatches    = 100;
  int test_type   = 3;
  int solver_type = 0;

  /* Parse command line arguments and setup UserData */
  int argi = 0;
  if (argc > 1) { nbatches = atoi(argv[++argi]); }
  if (argc > 2) { solver_type = atoi(argv[++argi]); }
  if (argc > 3) { test_type = atoi(argv[++argi]); }

  UserData udata{nbatches, batchSize, memhelper};

  /* Set the Reaction parameters according to test_type_type */
  for (int batchj = 0; batchj < udata.nbatches; ++batchj)
  {
    if (test_type == 1)
    {
      udata.u0[batchj] = SUN_RCONST(3.9);
      udata.v0[batchj] = SUN_RCONST(1.1);
      udata.w0[batchj] = SUN_RCONST(2.8);
      udata.a[batchj]  = SUN_RCONST(1.2);
      udata.b[batchj]  = SUN_RCONST(2.5);
      udata.ep[batchj] = SUN_RCONST(1.0e-5);
    }
    else if (test_type == 2)
    {
      udata.u0[batchj] = SUN_RCONST(3.0);
      udata.v0[batchj] = SUN_RCONST(3.0);
      udata.w0[batchj] = SUN_RCONST(3.5);
      udata.a[batchj]  = SUN_RCONST(0.5);
      udata.b[batchj]  = SUN_RCONST(3.0);
      udata.ep[batchj] = SUN_RCONST(5.0e-4);
    }
    else if (test_type == 3)
    {
      udata.u0[batchj] = SUN_RCONST(1.2);
      udata.v0[batchj] = SUN_RCONST(3.1);
      udata.w0[batchj] = SUN_RCONST(3.0);
      udata.a[batchj]  = SUN_RCONST(1.0);
      udata.b[batchj]  = SUN_RCONST(3.5);
      udata.ep[batchj] = SUN_RCONST(5.0e-6);
    }
  }

  /* Create CUDA or HIP vector for I.C. and abstol vector */
  y = HIP_OR_CUDA(N_VNew_Hip(batchSize * nbatches, sunctx);
                  , N_VNew_Cuda(batchSize * nbatches,
                                sunctx);) if (check_retval((void*)y, "N_VNew", 0))
  {
    return (1);
  }
  abstol = N_VClone(y);
  if (check_retval((void*)abstol, "N_VClone", 0)) { return (1); }

  /* Initialize y */
  sunrealtype* ydata = N_VGetArrayPointer(y);
  for (int batchj = 0; batchj < udata.nbatches; ++batchj)
  {
    ydata[batchj * udata.batchSize]     = udata.u0[batchj];
    ydata[batchj * udata.batchSize + 1] = udata.v0[batchj];
    ydata[batchj * udata.batchSize + 2] = udata.w0[batchj];
  }
  HIP_OR_CUDA(N_VCopyToDevice_Hip(y);, N_VCopyToDevice_Cuda(y);)

  /* Set the vector absolute tolerance */
  sunrealtype* abstol_data = N_VGetArrayPointer(abstol);
  for (int batchj = 0; batchj < udata.nbatches; ++batchj)
  {
    abstol_data[batchj * udata.batchSize]     = 1.0e-10;
    abstol_data[batchj * udata.batchSize + 1] = 1.0e-10;
    abstol_data[batchj * udata.batchSize + 2] = 1.0e-10;
  }
  HIP_OR_CUDA(N_VCopyToDevice_Hip(abstol);, N_VCopyToDevice_Cuda(abstol);)

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  /* Call CVodeSetUserData to attach the user data structure */
  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return (1); }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return (1); }

  /* Create the SUNLinearSolver object for use by CVode */
  if (solver_type == 0)
  {
    /* Create SUNMatrix for use in linear solves */
    A = SUNMatrix_MagmaDenseBlock(udata.nbatches, udata.batchSize,
                                  udata.batchSize, SUNMEMTYPE_DEVICE, memhelper,
                                  NULL, sunctx);
    if (check_retval((void*)A, "SUNMatrix_MagmaDenseBlock", 0)) { return (1); }

    LS = SUNLinSol_MagmaDense(y, A, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_MagmaDense", 0)) { return (1); }

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(cvode_mem, Jac);
    if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }
  }
  else
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_MagmaDense", 0)) { return (1); }

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }
  }

  /* In loop, call CVode, print results, and test_type for error.
     Break out of loop when preset output times have been reached.  */
  printf(" \nBatch of independent 3-species kinetics problems\n\n");
  printf("number of batches = %d\n\n", udata.nbatches);

  sunrealtype t    = T0;
  sunrealtype tout = T0 + dTout;
  for (int iout = 0; iout < Nt; iout++)
  {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    HIP_OR_CUDA(N_VCopyFromDevice_Hip(y);, N_VCopyFromDevice_Cuda(y);)
    for (int batchj = 0; batchj < udata.nbatches; batchj += 10)
    {
      printf("batch_%d: ", batchj);
      PrintOutput(t, ydata[batchj * udata.batchSize],
                  ydata[batchj * udata.batchSize + 1],
                  ydata[batchj * udata.batchSize + 2]);
    }

    if (check_retval(&retval, "CVode", 1)) { break; }
    if (retval == CV_SUCCESS)
    {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem, LS);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  if (A) { SUNMatDestroy(A); }

  return (0);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* Right hand side function. This just launches the CUDA or HIP kernel
   to do the actual computation. At the very least, doing this
   saves moving the vector data in y and ydot to/from the device
   every evaluation of f. */
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata;
  sunrealtype *ydata, *ydotdata;

  udata    = (UserData*)user_data;
  ydata    = N_VGetDeviceArrayPointer(y);
  ydotdata = N_VGetDeviceArrayPointer(ydot);

  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->nbatches + threads_per_block - 1) /
                        threads_per_block;
  f_kernel<<<num_blocks, threads_per_block>>>(t, ydata, ydotdata, udata->a.get(),
                                              udata->b.get(), udata->ep.get(),
                                              udata->nbatches, udata->batchSize);

  HIP_OR_CUDA(hipDeviceSynchronize();, cudaDeviceSynchronize();)
  HIP_OR_CUDA(hipError_t cuerr = hipGetLastError();
              , cudaError_t cuerr = cudaGetLastError();)
  if (cuerr != HIP_OR_CUDA(hipSuccess, cudaSuccess))
  {
    fprintf(stderr, ">>> ERROR in f: cudaGetLastError returned %s\n",
            HIP_OR_CUDA(hipGetErrorName(cuerr), cudaGetErrorName(cuerr)));
    return (-1);
  }

  return (0);
}

/* Right hand side function evaluation kernel.
   This kernel needs threads >= nbatches. */
__global__ void f_kernel(sunrealtype t, sunrealtype* ydata,
                         sunrealtype* ydotdata, sunrealtype* A, sunrealtype* B,
                         sunrealtype* Ep, int nbatches, int batchSize)
{
  sunrealtype u, v, w, a, b, ep;

  for (int batchj = blockIdx.x * blockDim.x + threadIdx.x; batchj < nbatches;
       batchj += blockDim.x * gridDim.x)
  {
    a = A[batchj];
    b = B[batchj], ep = Ep[batchj];

    u = ydata[batchj * batchSize];
    v = ydata[batchj * batchSize + 1];
    w = ydata[batchj * batchSize + 2];

    ydotdata[batchj * batchSize]     = a - (w + 1.0) * u + v * u * u;
    ydotdata[batchj * batchSize + 1] = w * u - v * u * u;
    ydotdata[batchj * batchSize + 2] = (b - w) / ep - w * u;
  }
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
 * This is done on the GPU.
 */

int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata = (UserData*)user_data;

  sunrealtype* Jdata = SUNMatrix_MagmaDense_Data(J);
  sunrealtype* ydata = N_VGetDeviceArrayPointer(y);

  unsigned block_size = HIP_OR_CUDA(64, 32);
  unsigned grid_size  = (udata->nbatches + block_size - 1) / block_size;

  j_kernel<<<grid_size, block_size>>>(ydata, Jdata, udata->a.get(),
                                      udata->b.get(), udata->ep.get(),
                                      udata->nbatches, udata->batchSize);

  HIP_OR_CUDA(hipDeviceSynchronize();, cudaDeviceSynchronize();)
  HIP_OR_CUDA(hipError_t cuerr = hipGetLastError();
              , cudaError_t cuerr = cudaGetLastError();)
  if (cuerr != HIP_OR_CUDA(hipSuccess, cudaSuccess))
  {
    fprintf(stderr, ">>> ERROR in Jac: cudaGetLastError returned %s\n",
            HIP_OR_CUDA(hipGetErrorName(cuerr), cudaGetErrorName(cuerr)));
    return (-1);
  }

  return (0);
}

/* Jacobian evaluation GPU kernel */
__global__ void j_kernel(sunrealtype* ydata, sunrealtype* Jdata, sunrealtype* A,
                         sunrealtype* B, sunrealtype* Ep, int nbatches,
                         int batchSize)
{
  int N  = batchSize;
  int NN = N * N;
  sunrealtype u, v, w, ep;

  for (int batchj = blockIdx.x * blockDim.x + threadIdx.x; batchj < nbatches;
       batchj += blockDim.x * gridDim.x)
  {
    ep = Ep[batchj];

    /* get y values */
    u = ydata[N * batchj];
    v = ydata[N * batchj + 1];
    w = ydata[N * batchj + 2];

    /* first col of block */
    Jdata[NN * batchj]     = -(w + 1.0) + 2.0 * u * v;
    Jdata[NN * batchj + 1] = u * u;
    Jdata[NN * batchj + 2] = -u;

    /* second col of block */
    Jdata[NN * batchj + 3] = u * u;
    Jdata[NN * batchj + 4] = -u * u;
    Jdata[NN * batchj + 5] = u;

    /* third col of block */
    Jdata[NN * batchj + 6] = -w;
    Jdata[NN * batchj + 7] = 0.0;
    Jdata[NN * batchj + 8] = -1.0 / ep - u;
  }
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2, sunrealtype y3)
{
#if defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

/*
 * Get and print some final statistics
 */

void PrintFinalStats(void* cvode_mem, SUNLinearSolver LS)
{
  long int nst, nfe, nsetups, nje, nni, ncfn, netf, nge;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n", nst, nfe,
         nsetups, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n \n", nni, ncfn,
         netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s(failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}
