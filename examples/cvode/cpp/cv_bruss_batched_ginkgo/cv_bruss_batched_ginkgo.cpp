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
 * We simulate a scenario where a set of independent ODEs are batched together
 * to form a larger system. Each independent ODE system has 3 components,
 * Y = [u,v,w], satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions Y0 = [u0,v0,w0].
 * We replicate this ODE system num_batches times, with each batch using a
 * different value of ep, to eventually construct an ODE system with
 * 3*num_batches components.
 *
 * The problem is stiff. We have 3 different reactors for testing:
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
 * systems results in a batch diagonal linear system, with the Ginkgo
 * SUNLinearSolver which supports batched iterative methods. 100 outputs are
 * printed at equal intervals, and run statistics are printed at the end.
 *
 * The program takes three optional arguments, the number of independent ODE
 * systems (i.e. number of batches), the linear solver type
 * (ginkgo batched BiCGSTAB, non-batched SUNDIALS BiCGSTAB with the Jaocbian
 * computed by difference quotients, or non-batched BiCGSTAB with an analytical Jacobian)
 * , and the test type (uniform_1, uniform_2, uniform_3, or random).
 *
 *    ./cv_bruss_batched_ginkgo [num_batches] [solver_type] [test_type]
 *
 * Options:
 *    num_batches <int>
 *    solver_type:
 *       0 - Ginkgo batched BiCGSTAB
 *       1 - SUNDIALS non-batched BiCGSTAB with Jacobian-vector products computed by difference quotients
 *       2 - SUNDIALS non-batched BiCGSTAB with analytical Jacobian-vector products
 *    test_type:
 *       0 - all batches are reactor type 0
 *       1 - all batches are reactor type 1
 *       2 - all batches are reactor type 2
 *       3 - random distribution of reactors
 *       4 - all batches are reactor type 0 but with random distribution of ep (stiffness)
 *       5 - all batches are reactor type 1 but with random distribution of ep (stiffness)
 *       6 - all batches are reactor type 2 but with random distribution of ep (stiffness)
 *       7 - random distribution of reactors and ep (stiffness)
 * --------------------------------------------------------------------------*/

#include <cstdio>
#include <memory>
#include <random>

#include <cvode/cvode.h> /* prototypes for CVODE fcts., consts. */
#include <sunlinsol/sunlinsol_ginkgobatch.hpp> /* access to Ginkgo batch SUNLinearSolver */
#include <sunlinsol/sunlinsol_spbcgs.h> /* access to SUNDIALS BiCGSTAB SUNLinearSolver */

#include "cv_bruss_batched_ginkgo.hpp"

/* Functions called by the integrator */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int JacVec(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* user_data, N_Vector tmp);

/* Private functions which implement the actual operators either on the CPU or GPU */
SERIAL_CUDA_OR_HIP(, __global__, __global__)
static void f_kernel(sunrealtype t, sunrealtype* y, sunrealtype* ydot,
                     sunrealtype* A, sunrealtype* B, sunrealtype* Ep, int neq,
                     int num_batches, int batch_size);
SERIAL_CUDA_OR_HIP(, __global__, __global__)
static void j_kernel(sunrealtype* ydata, sunrealtype* Jdata, sunrealtype* A,
                     sunrealtype* B, sunrealtype* Ep, int neq, int num_batches,
                     int batch_size, int nnzper);

SERIAL_CUDA_OR_HIP(, __global__, __global__)
static void jv_kernel(sunrealtype* vdata, sunrealtype* Jvdata, sunrealtype* ydata,
                      sunrealtype* A, sunrealtype* B, sunrealtype* Ep, int neq,
                      int num_batches, int batch_size, int nnzper);

/* Private function to initialize the Jacobian sparsity pattern */
static int JacInit(SUNMatrix J);

/* Private function to output results */
static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3);

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Shortcuts */
using GkoBatchMatrixType = gko::batch::matrix::Csr<sunrealtype, sunindextype>;
using GkoBatchSolverType = gko::batch::solver::Bicgstab<sunrealtype>;
using SUNGkoMatrixType   = sundials::ginkgo::BatchMatrix<GkoBatchMatrixType>;
using SUNGkoLinearSolverType =
  sundials::ginkgo::BatchLinearSolver<GkoBatchSolverType, GkoBatchMatrixType>;

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char* argv[])
{
  const sunrealtype T0 = SUN_RCONST(0.0);    /* initial time                  */
  const sunrealtype Tf = SUN_RCONST(10.0);   /* final time                    */
  const sunrealtype dTout = SUN_RCONST(5.0); /* time between outputs          */
  const int Nt = (int)ceil(Tf / dTout);      /* number of output times        */
  const sunrealtype reltol = SUN_RCONST(1.0e-6); /* relative integrator tolerance */
  int retval;

  /* Create the SUNDIALS context */
  sundials::Context sunctx;

  /* Create the Ginkgo executor */
#if defined(USE_CUDA)
  auto gko_exec{gko::CudaExecutor::create(0, gko::OmpExecutor::create())};
#elif defined(USE_HIP)
  auto gko_exec{gko::HipExecutor::create(0, gko::OmpExecutor::create())};
#elif defined(USE_OMP)
  auto gko_exec{gko::OmpExecutor::create()};
#else
  auto gko_exec{gko::ReferenceExecutor::create()};
#endif

  /* Create a SUNMemoryHelper for HIP or CUDA depending on the target.
     This will be used underneath the Arrays. */
  SUNMemoryHelper memhelper = SERIAL_CUDA_OR_HIP(SUNMemoryHelper_Sys(sunctx),
                                                 SUNMemoryHelper_Cuda(sunctx),
                                                 SUNMemoryHelper_Hip(sunctx));

  /* Set defaults */
  int batch_size  = 3;
  int num_batches = 100;
  int test_type   = 0;
  int solver_type = 0;

  /* Parse command line arguments and setup UserData */
  int argi = 0;
  if (argc > 1) { num_batches = atoi(argv[++argi]); }
  if (argc > 2) { solver_type = atoi(argv[++argi]); }
  if (argc > 3) { test_type = atoi(argv[++argi]); }

  UserData udata{num_batches, batch_size, batch_size * batch_size, memhelper};

  /* Set the Reaction parameters according to test_type */
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<sunrealtype> dist(SUN_RCONST(0.0),
                                                   SUN_RCONST(8.0));
  std::uniform_int_distribution<> dist_int(0, 2);

  bool random_reactors  = test_type == 3 || test_type == 7;
  bool random_stiffness = test_type == 4 || test_type == 5 || test_type == 6 ||
                          test_type == 7;
  for (int batchj = 0; batchj < udata.num_batches; ++batchj)
  {
    int reactor_type = random_reactors ? dist_int(gen) : test_type % 4;
    if (reactor_type == 0)
    {
      udata.u0[batchj] = SUN_RCONST(3.9);
      udata.v0[batchj] = SUN_RCONST(1.1);
      udata.w0[batchj] = SUN_RCONST(2.8);
      udata.a[batchj]  = SUN_RCONST(1.2);
      udata.b[batchj]  = SUN_RCONST(2.5);
      udata.ep[batchj] = random_stiffness
                           ? std::pow(SUN_RCONST(10.0), -dist(gen))
                           : SUN_RCONST(1.0e-5);
    }
    else if (reactor_type == 1)
    {
      udata.u0[batchj] = SUN_RCONST(3.0);
      udata.v0[batchj] = SUN_RCONST(3.0);
      udata.w0[batchj] = SUN_RCONST(3.5);
      udata.a[batchj]  = SUN_RCONST(0.5);
      udata.b[batchj]  = SUN_RCONST(3.0);
      udata.ep[batchj] = random_stiffness
                           ? std::pow(SUN_RCONST(10.0), -dist(gen))
                           : SUN_RCONST(5.0e-4);
    }
    else if (reactor_type == 2)
    {
      udata.u0[batchj] = SUN_RCONST(1.2);
      udata.v0[batchj] = SUN_RCONST(3.1);
      udata.w0[batchj] = SUN_RCONST(3.0);
      udata.a[batchj]  = SUN_RCONST(1.0);
      udata.b[batchj]  = SUN_RCONST(3.5);
      udata.ep[batchj] = random_stiffness
                           ? std::pow(SUN_RCONST(10.0), -dist(gen))
                           : SUN_RCONST(5.0e-6);
    }
  }

  /* Create CUDA or HIP vector of length neq for I.C. and abstol vector */
  N_Vector y = N_VNew(udata.neq, sunctx);
  if (check_retval((void*)y, "N_VNew", 0)) { return 1; }

  N_Vector abstol = N_VClone(y);
  if (check_retval((void*)abstol, "N_VClone", 0)) { return 1; }

  /* Initialize y */
  sunrealtype* ydata = N_VGetArrayPointer(y);
  for (int batchj = 0; batchj < udata.num_batches; ++batchj)
  {
    ydata[batchj * udata.batch_size]     = udata.u0[batchj];
    ydata[batchj * udata.batch_size + 1] = udata.v0[batchj];
    ydata[batchj * udata.batch_size + 2] = udata.w0[batchj];
  }

  /* Set the vector absolute tolerance */
  sunrealtype* abstol_data = N_VGetArrayPointer(abstol);
  for (int batchj = 0; batchj < udata.num_batches; ++batchj)
  {
    abstol_data[batchj * udata.batch_size]     = 1.0e-15;
    abstol_data[batchj * udata.batch_size + 1] = 1.0e-15;
    abstol_data[batchj * udata.batch_size + 2] = 1.0e-15;
  }

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return 1; }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return 1; }

  /* Call CVodeSetUserData to attach the user data structure */
  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return 1; }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return 1; }

  /* Create SUNMatrix for use in linear solves */
  sunindextype batch_rows = udata.batch_size;
  sunindextype batch_cols = udata.batch_size;
  auto batch_dim          = gko::batch_dim<2>(udata.num_batches,
                                     gko::dim<2>(batch_rows, batch_cols));
  auto gko_matA =
    gko::share(GkoBatchMatrixType::create(gko_exec, batch_dim, udata.nnzper));
  SUNGkoMatrixType A{gko_matA, sunctx};

  /* Initialize the Jacobian with its fixed sparsity pattern */
  retval = JacInit(A);
  if (check_retval(&retval, "JacInit", 1)) { return 1; }

  /* Create the SUNLinearSolver object for use by CVode */
  SUNGkoLinearSolverType LS{gko_exec,
                            gko::batch::stop::tolerance_type::absolute, nullptr,
                            static_cast<gko::size_type>(udata.num_batches),
                            sunctx};

  SUNLinearSolver SUNLS = NULL;

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  if (solver_type == 0)
  {
    retval = CVodeSetLinearSolver(cvode_mem, LS.get(), A.get());
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return 1; }

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(cvode_mem, Jac);
    if (check_retval(&retval, "CVodeSetJacFn", 1)) return 1;

    /* WARNING: THIS STEP IS ESSENTIAL WHEN USING SUNLINSOL_GINKGOBATCH!
     We need to adjust the norm factor to be the sqrt of the length of the batch
     instead of the sqrt of the length of the overall system, since the linear
     solver is batched.
    */
    retval = CVodeSetLSNormFactor(cvode_mem, std::sqrt(udata.batch_size));
    if (check_retval(&retval, "CVodeSetJacFn", 1)) return 1;
  }
  else
  {
    SUNLS  = SUNLinSol_SPBCGS(y, SUN_PREC_NONE, 0, sunctx);
    retval = CVodeSetLinearSolver(cvode_mem, SUNLS, NULL);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return 1; }

    if (solver_type == 2)
    {
      /* Set the user-supplied Jacobian-vector product routine */
      retval = CVodeSetJacTimes(cvode_mem, NULL, JacVec);
      if (check_retval(&retval, "CVodeSetJacTimes", 1)) return 1;
    }
  }

  /* Increase the max number of time steps allowed */
  retval = CVodeSetMaxNumSteps(cvode_mem, 5000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return 1;

  /* In loop, call CVode, print results, and test_type for error.
     Break out of loop when preset output times have been reached.  */
  printf(" \nBatch of independent 3-species kinetics problems\n");
  printf("number of batches = %d, solver_type = %d, test type = %d\n\n",
         udata.num_batches, solver_type, test_type);

  sunrealtype t    = T0;
  sunrealtype tout = T0 + dTout;
  for (int iout = 0; iout < Nt; iout++)
  {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    for (int batchj = 0; batchj < udata.num_batches; batchj++)
    {
      printf("batch_%d: ", batchj);
      PrintOutput(t, ydata[batchj * udata.batch_size],
                  ydata[batchj * udata.batch_size + 1],
                  ydata[batchj * udata.batch_size + 2]);
    }

    if (check_retval(&retval, "CVode", 1)) break;
    if (retval == CV_SUCCESS)
    {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }

  /* Print some final statistics */
  printf("Final integrator statistics:\n");
  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "CVodePrintAllStats", 1)) return 1;

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free linear solver if needed */
  if (SUNLS) SUNLinSolFree(SUNLS);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return 0;
}

/*
 * Jacobian initialization routine. This sets the sparisty pattern of
 * the batches of the Jacobian J(t,y) = df/dy. This is performed on the CPU,
 * and only occurs at the beginning of the simulation.
 */

int JacInit(SUNMatrix J)
{
  auto Jgko     = static_cast<SUNGkoMatrixType*>(J->content)->GkoMtx();
  auto gko_exec = Jgko->get_executor();

  sunindextype* rowptrs = Jgko->get_row_ptrs();
  if (rowptrs == nullptr) { return -1; }
  sunindextype* colvals = Jgko->get_col_idxs();
  if (colvals == nullptr) { return -1; }

#if defined(USE_CUDA) || defined(USE_HIP)

  sunindextype rowptrs_host[4];

  /* there are 3 entries per row */
  rowptrs_host[0] = 0;
  rowptrs_host[1] = 3;
  rowptrs_host[2] = 6;
  rowptrs_host[3] = 9;

  gko_exec->copy_from(gko_exec->get_master(), 4, rowptrs_host, rowptrs);

  sunindextype colvals_host[9];

  /* first row of batch */
  colvals_host[0] = 0;
  colvals_host[1] = 1;
  colvals_host[2] = 2;

  /* second row of batch */
  colvals_host[3] = 0;
  colvals_host[4] = 1;
  colvals_host[5] = 2;

  /* third row of batch */
  colvals_host[6] = 0;
  colvals_host[7] = 1;
  colvals_host[8] = 2;

  gko_exec->copy_from(gko_exec->get_master(), 9, colvals_host, colvals);

#else

  /* there are 3 entries per row */
  rowptrs[0] = 0;
  rowptrs[1] = 3;
  rowptrs[2] = 6;
  rowptrs[3] = 9;

  /* first row of batch */
  colvals[0] = 0;
  colvals[1] = 1;
  colvals[2] = 2;

  /* second row of batch */
  colvals[3] = 0;
  colvals[4] = 1;
  colvals[5] = 2;

  /* third row of batch */
  colvals[6] = 0;
  colvals[7] = 1;
  colvals[8] = 2;

#endif

  return 0;
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
  ydata    = N_VGetArrayPointer(y);
  ydotdata = N_VGetArrayPointer(ydot);

#if defined(USE_CUDA) || defined(USE_HIP)
  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->num_batches + threads_per_block - 1) /
                        threads_per_block;
  f_kernel<<<num_blocks, threads_per_block>>>(t, ydata, ydotdata, udata->a.get(),
                                              udata->b.get(), udata->ep.get(),
                                              udata->neq, udata->num_batches,
                                              udata->batch_size);

  SERIAL_CUDA_OR_HIP(, cudaDeviceSynchronize(), hipDeviceSynchronize());
  SERIAL_CUDA_OR_HIP(, cudaError_t gpu_err = cudaGetLastError(),
                     hipError_t gpu_err    = hipGetLastError());
  if (gpu_err != SERIAL_CUDA_OR_HIP(, cudaSuccess, hipSuccess))
  {
    fprintf(stderr, ">>> ERROR in f: GetLastError returned %s\n",
            SERIAL_CUDA_OR_HIP(, cudaGetErrorName(gpu_err),
                               hipGetErrorName(gpu_err)));
    return -1;
  }
#else
  f_kernel(t, ydata, ydotdata, udata->a.get(), udata->b.get(), udata->ep.get(),
           udata->neq, udata->num_batches, udata->batch_size);
#endif

  return 0;
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
 * This is done on the GPU.
 */

int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata = (UserData*)user_data;
  auto Jgko       = static_cast<SUNGkoMatrixType*>(J->content)->GkoMtx();

  sunrealtype* Jdata = Jgko->get_values();
  sunrealtype* ydata = N_VGetArrayPointer(y);

#if defined(USE_CUDA) || defined(USE_HIP)
  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->num_batches + threads_per_block - 1) /
                        threads_per_block;

  j_kernel<<<num_blocks, threads_per_block>>>(ydata, Jdata, udata->a.get(),
                                              udata->b.get(), udata->ep.get(),
                                              udata->neq, udata->num_batches,
                                              udata->batch_size, udata->nnzper);

  SERIAL_CUDA_OR_HIP(, cudaDeviceSynchronize(), hipDeviceSynchronize());
  SERIAL_CUDA_OR_HIP(, cudaError_t gpu_err = cudaGetLastError(),
                     hipError_t gpu_err    = hipGetLastError());
  if (gpu_err != SERIAL_CUDA_OR_HIP(, cudaSuccess, hipSuccess))
  {
    fprintf(stderr, ">>> ERROR in Jac: GetLastError returned %s\n",
            SERIAL_CUDA_OR_HIP(, cudaGetErrorName(gpu_err),
                               hipGetErrorName(gpu_err)));
    return -1;
  }
#else
  j_kernel(ydata, Jdata, udata->a.get(), udata->b.get(), udata->ep.get(),
           udata->neq, udata->num_batches, udata->batch_size, udata->nnzper);
#endif

  return 0;
}

int JacVec(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy,
           void* user_data, N_Vector tmp)
{
  UserData* udata = (UserData*)user_data;

  sunrealtype* vdata  = N_VGetArrayPointer(v);
  sunrealtype* Jvdata = N_VGetArrayPointer(Jv);
  sunrealtype* ydata  = N_VGetArrayPointer(y);

#if defined(USE_CUDA) || defined(USE_HIP)
  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->num_batches + threads_per_block - 1) /
                        threads_per_block;

  jv_kernel<<<num_blocks, threads_per_block>>>(vdata, Jvdata, ydata,
                                               udata->a.get(), udata->b.get(),
                                               udata->ep.get(), udata->neq,
                                               udata->num_batches,
                                               udata->batch_size, udata->nnzper);

  SERIAL_CUDA_OR_HIP(, cudaDeviceSynchronize(), hipDeviceSynchronize());
  SERIAL_CUDA_OR_HIP(, cudaError_t gpu_err = cudaGetLastError(),
                     hipError_t gpu_err    = hipGetLastError());
  if (gpu_err != SERIAL_CUDA_OR_HIP(, cudaSuccess, hipSuccess))
  {
    fprintf(stderr, ">>> ERROR in JacVec: GetLastError returned %s\n",
            SERIAL_CUDA_OR_HIP(, cudaGetErrorName(gpu_err),
                               hipGetErrorName(gpu_err)));
    return -1;
  }
#else
  jv_kernel(vdata, Jvdata, ydata, udata->a.get(), udata->b.get(), udata->ep.get(),
            udata->neq, udata->num_batches, udata->batch_size, udata->nnzper);
#endif

  return 0;
}

#if defined(USE_CUDA) || defined(USE_HIP)
/* Right hand side function evaluation GPU kernel.
   This kernel needs total number of threads >= num_batches. */
__global__ void f_kernel(sunrealtype t, sunrealtype* ydata, sunrealtype* ydotdata,
                         sunrealtype* A, sunrealtype* B, sunrealtype* Ep,
                         int neq, int num_batches, int batch_size)
{
  sunrealtype u, v, w, a, b, ep;

  int batchj = blockIdx.x * blockDim.x + threadIdx.x;

  if (batchj < num_batches)
  {
    a = A[batchj];
    b = B[batchj], ep = Ep[batchj];

    u = ydata[batchj * batch_size];
    v = ydata[batchj * batch_size + 1];
    w = ydata[batchj * batch_size + 2];

    ydotdata[batchj * batch_size]     = a - (w + 1.0) * u + v * u * u;
    ydotdata[batchj * batch_size + 1] = w * u - v * u * u;
    ydotdata[batchj * batch_size + 2] = (b - w) / ep - w * u;
  }
}

/* Jacobian evaluation GPU kernel
   This kernel needs total number of threads >= num_batches. */
__global__ void j_kernel(sunrealtype* ydata, sunrealtype* Jdata, sunrealtype* A,
                         sunrealtype* B, sunrealtype* Ep, int neq,
                         int num_batches, int batch_size, int nnzper)
{
  sunrealtype u, v, w, ep;

  int batchj = blockIdx.x * blockDim.x + threadIdx.x;

  if (batchj < num_batches)
  {
    ep = Ep[batchj];

    /* get y values */
    u = ydata[batch_size * batchj];
    v = ydata[batch_size * batchj + 1];
    w = ydata[batch_size * batchj + 2];

    /* first row of batch */
    Jdata[nnzper * batchj]     = -(w + 1.0) + 2.0 * u * v;
    Jdata[nnzper * batchj + 1] = u * u;
    Jdata[nnzper * batchj + 2] = -u;

    /* second row of batch */
    Jdata[nnzper * batchj + 3] = w - 2.0 * u * v;
    Jdata[nnzper * batchj + 4] = -u * u;
    Jdata[nnzper * batchj + 5] = u;

    /* third row of batch */
    Jdata[nnzper * batchj + 6] = -w;
    Jdata[nnzper * batchj + 7] = 0.0;
    Jdata[nnzper * batchj + 8] = -1.0 / ep - u;
  }
}

/* Jacobian-vector product GPU kernel
   This kernel needs total number of threads >= num_batches. */
__global__ void jv_kernel(sunrealtype* vdata, sunrealtype* Jvdata,
                          sunrealtype* ydata, sunrealtype* A, sunrealtype* B,
                          sunrealtype* Ep, int neq, int num_batches,
                          int batch_size, int nnzper)
{
  sunrealtype u, v, w, v0, v1, v2, ep;

  int batchj = blockIdx.x * blockDim.x + threadIdx.x;

  if (batchj < num_batches)
  {
    ep = Ep[batchj];

    /* get y values */
    u = ydata[batch_size * batchj];
    v = ydata[batch_size * batchj + 1];
    w = ydata[batch_size * batchj + 2];

    /* get v values */
    v0 = vdata[batch_size * batchj];
    v1 = vdata[batch_size * batchj + 1];
    v2 = vdata[batch_size * batchj + 2];

    /* initialize Jv to zero */
    Jvdata[batch_size * batchj]     = 0.0;
    Jvdata[batch_size * batchj + 1] = 0.0;
    Jvdata[batch_size * batchj + 2] = 0.0;

    // add J[:,0]*v[0]
    Jvdata[batch_size * batchj] += (-(w + 1.0) + 2.0 * u * v) * v0;
    Jvdata[batch_size * batchj + 1] += (w - 2.0 * u * v) * v0;
    Jvdata[batch_size * batchj + 2] += -w * v0;

    // add J[:,1]*v[1]
    Jvdata[batch_size * batchj] += (u * u) * v1;
    Jvdata[batch_size * batchj + 1] += (-u * u) * v1;
    Jvdata[batch_size * batchj + 2] += 0.0;

    // add J[:,2]*v[2]
    Jvdata[batch_size * batchj] += -u * v2;
    Jvdata[batch_size * batchj + 1] += u * v2;
    Jvdata[batch_size * batchj + 2] += (-1.0 / ep - u) * v2;
  }
}
#else
void f_kernel(sunrealtype t, sunrealtype* ydata, sunrealtype* ydotdata,
              sunrealtype* A, sunrealtype* B, sunrealtype* Ep, int neq,
              int num_batches, int batch_size)
{
  sunrealtype u, v, w, a, b, ep;

  for (int batchj = 0; batchj < num_batches; batchj++)
  {
    a = A[batchj];
    b = B[batchj], ep = Ep[batchj];

    u = ydata[batchj * batch_size];
    v = ydata[batchj * batch_size + 1];
    w = ydata[batchj * batch_size + 2];

    ydotdata[batchj * batch_size]     = a - (w + 1.0) * u + v * u * u;
    ydotdata[batchj * batch_size + 1] = w * u - v * u * u;
    ydotdata[batchj * batch_size + 2] = (b - w) / ep - w * u;
  }
}

void j_kernel(sunrealtype* ydata, sunrealtype* Jdata, sunrealtype* A,
              sunrealtype* B, sunrealtype* Ep, int neq, int num_batches,
              int batch_size, int nnzper)
{
  sunrealtype u, v, w, ep;
  for (int batchj = 0; batchj < num_batches; batchj++)
  {
    ep = Ep[batchj];
    u  = ydata[batch_size * batchj];
    v  = ydata[batch_size * batchj + 1];
    w  = ydata[batch_size * batchj + 2];
    // first row of batch
    Jdata[nnzper * batchj]     = -(w + 1.0) + 2.0 * u * v;
    Jdata[nnzper * batchj + 1] = u * u;
    Jdata[nnzper * batchj + 2] = -u;
    // second row of batch
    Jdata[nnzper * batchj + 3] = w - 2.0 * u * v;
    Jdata[nnzper * batchj + 4] = -u * u;
    Jdata[nnzper * batchj + 5] = u;
    // third row of batch
    Jdata[nnzper * batchj + 6] = -w;
    Jdata[nnzper * batchj + 7] = 0.0;
    Jdata[nnzper * batchj + 8] = -1.0 / ep - u;
  }
}

void jv_kernel(sunrealtype* vdata, sunrealtype* Jvdata, sunrealtype* ydata,
               sunrealtype* A, sunrealtype* B, sunrealtype* Ep, int neq,
               int num_batches, int batch_size, int nnzper)
{
  sunrealtype u, v, w, v0, v1, v2, ep;
  for (int batchj = 0; batchj < num_batches; batchj++)
  {
    ep = Ep[batchj];

    // get y values
    u = ydata[batch_size * batchj];
    v = ydata[batch_size * batchj + 1];
    w = ydata[batch_size * batchj + 2];

    // get v values
    v0 = vdata[batch_size * batchj];
    v1 = vdata[batch_size * batchj + 1];
    v2 = vdata[batch_size * batchj + 2];

    // initialize Jv to zero
    Jvdata[batch_size * batchj]     = 0.0;
    Jvdata[batch_size * batchj + 1] = 0.0;
    Jvdata[batch_size * batchj + 2] = 0.0;

    // add J[:,0]*v[0]
    Jvdata[batch_size * batchj] += (-(w + 1.0) + 2.0 * u * v) * v0;
    Jvdata[batch_size * batchj + 1] += (w - 2.0 * u * v) * v0;
    Jvdata[batch_size * batchj + 2] += -w * v0;

    // add J[:,1]*v[1]
    Jvdata[batch_size * batchj] += (u * u) * v1;
    Jvdata[batch_size * batchj + 1] += (-u * u) * v1;
    Jvdata[batch_size * batchj + 2] += 0.0;

    // add J[:,2]*v[2]
    Jvdata[batch_size * batchj] += -u * v2;
    Jvdata[batch_size * batchj + 1] += u * v2;
    Jvdata[batch_size * batchj + 2] += (-1.0 / ep - u) * v2;
  }
}
#endif

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2, sunrealtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
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
