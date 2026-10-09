/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                David J. Gardner and Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * The following is a simple example problem based off of ark_brusselator.c.
 *
 * We simulate a scenario where a set of independent ODEs are batched together
 * to form a larger system. Each independent ODE system has 3 components,
 * Y = [u, v, w], satisfying the equations,
 *
 *   du/dt = a - (w + 1) * u + v * u^2
 *   dv/dt = w * u - v * u^2
 *   dw/dt = (b - w) / ep - w * u
 *
 * for t in the interval [0, 10], with initial conditions Y0 = [u0, v0, w0].
 * The problem is stiff and there are 3 testing scenarios:
 *
 * Reactor 0: u0 = 3.9, v0 = 1.1, w0 = 2.8, a = 1.2, b = 2.5, ep = 1.0e-5
 *   Here, all three components exhibit a rapid transient change during the
 *   first 0.2 time units, followed by a slow and smooth evolution.
 *
 * Reactor 1: u0 = 3, v0 = 3, w0 = 3.5, a = 0.5, b = 3, ep = 5.0e-4
 *   Here, all components undergo very rapid initial transients during the first
 *   0.3 time units, and all then proceed very smoothly for the remainder of the
 *   simulation.
 *
 * Reactor 2: u0 = 1.2, v0 = 3.1, w0 = 3, a = 1, b = 3.5, ep = 5.0e-6
 *   Here, w experiences a fast initial transient, jumping 0.5 within a few
 *   steps. All values proceed smoothly until around t=6.5, when both u and v
 *   undergo a sharp transition, with u increasing from around 0.5 to 5 and v
 *   decreasing from around 6 to 1 in less than 0.5 time units. After this
 *   transition, both u and v continue to evolve somewhat rapidly for another
 *   1.4 time units, and finish off smoothly.
 *
 * This program solves the problem with the BDF method, Newton iteration, a
 * user-supplied Jacobian routine, and, since the grouping of the independent
 * systems results in a block diagonal linear system, the dense KOKKOS
 * SUNLinearSolver which supports batched systems. 100 outputs are printed at
 * equal intervals, and run statistics are printed at the end.
 *
 * Unlike the example cv_bruss_batched_kokkos.cpp, this example utilizes Kokkos'
 * multi-dimensional view functionality to consider a 2D grouping, y(i,j), where
 * i corresponds with the batch index, and j corresponds to the component (u,v,w).
 *
 * The program takes three optional arguments, the number of independent ODE
 * systems (i.e., number of batches), the linear solver type (KOKKOS batched LU
 * or non-batched GMRES with the Jacobian computed by difference quotients)
 * the test type (uniform_0, uniform_1, or  uniform_2).
 *
 *   ./cv_bruss_batched_kokkos [num_batches] [solver_type] [test_type]
 *
 * Options:
 *   num_batches <int>
 *   solver_type:
 *     0 - KOKKOS batched LU (default)
 *     1 - SUNDIALS non-batched GMRES with difference quotients Jacobian
 *   test_type:
 *     0 - uniform_0, all batches are Reactor 0
 *     1 - uniform 1, all batches are Reactor 1
 *     2 - uniform 2, all batches are Reactor 2 (default)
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cvode/cvode.h>
#include <memory>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_kokkosdense.hpp>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>
#include <vector>

// Common utility functions
#include <example_utilities.hpp>

// Execution space
#if defined(USE_CUDA)
using ExecSpace = Kokkos::Cuda;
using MemSpace  = Kokkos::CudaSpace;
#elif defined(USE_HIP)
#if KOKKOS_VERSION / 10000 > 3
using ExecSpace = Kokkos::HIP;
using MemSpace  = Kokkos::HIPSpace;
#else
using ExecSpace = Kokkos::Experimental::HIP;
using MemSpace  = Kokkos::Experimental::HIPSpace;
#endif
#elif defined(USE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
using MemSpace  = Kokkos::HostSpace;
#else
using ExecSpace = Kokkos::Serial;
using MemSpace  = Kokkos::HostSpace;
#endif

using Vec1D     = Kokkos::View<sunrealtype*, MemSpace>;
using Vec2D     = Kokkos::View<sunrealtype**, Kokkos::LayoutRight, MemSpace>;
using Vec2DHost = Vec2D::HostMirror;
using VecType   = sundials::kokkos::Vector<ExecSpace>;
using MatType   = sundials::kokkos::DenseMatrix<ExecSpace>;
using LSType    = sundials::kokkos::DenseLinearSolver<ExecSpace>;
using SizeType  = VecType::size_type;

// Constants
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

// User-supplied functions called by CVODE
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// User data structure available in user-supplied callback functions
struct UserData
{
  int nbatches  = 100; // number of chemical networks
  int batchSize = 3;   // size of each network
  sunrealtype a, b;    // chemical concentrations that are constant
  sunrealtype ep;      // stiffness parameter
};

/* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
  // Create the SUNDIALS context
  sundials::Context sunctx;

  Kokkos::initialize(argc, argv);
  {
    // Create UserData
    UserData udata;

    // Parse command line options
    int argi = 0;

    // Total number of batch systems
    if (argc > 1) { udata.nbatches = atoi(argv[++argi]); }

    // Linear solver type
    int solver_type = 0;
    if (argc > 2) { solver_type = atoi(argv[++argi]); }

    // Problem setup
    int test_type = 2;
    if (argc > 3) { test_type = atoi(argv[++argi]); }

    // Shortcuts
    int nbatches  = udata.nbatches;
    int batchSize = udata.batchSize;

    std::cout << "\nBatch of independent 3-species kinetics problems\n"
              << "  number of batches = " << nbatches << "\n"
              << "  linear solver     = "
              << (solver_type ? "GMRES" : "KokkosKernels") << "\n"
              << "  test type         = " << test_type << "\n"
              << "  execution space   = " << ExecSpace().name() << "\n\n";

    sunrealtype u0, v0, w0;
    if (test_type == 0)
    {
      u0 = SUN_RCONST(3.9);
      v0 = SUN_RCONST(1.1);
      w0 = SUN_RCONST(2.8);

      udata.a  = SUN_RCONST(1.2);
      udata.b  = SUN_RCONST(2.5);
      udata.ep = SUN_RCONST(1.0e-5);
    }
    else if (test_type == 1)
    {
      u0 = SUN_RCONST(3.0);
      v0 = SUN_RCONST(3.0);
      w0 = SUN_RCONST(3.5);

      udata.a  = SUN_RCONST(0.5);
      udata.b  = SUN_RCONST(3.0);
      udata.ep = SUN_RCONST(5.0e-4);
    }
    else if (test_type == 2)
    {
      u0 = SUN_RCONST(1.2);
      v0 = SUN_RCONST(3.1);
      w0 = SUN_RCONST(3.0);

      udata.a  = SUN_RCONST(1.0);
      udata.b  = SUN_RCONST(3.5);
      udata.ep = SUN_RCONST(5.0e-6);
    }
    else
    {
      std::cerr << "ERROR: Invalid test type option\n";
      return -1;
    }

    // Create vector with the initial condition
    const sunrealtype T0 = SUN_RCONST(0.0);

    SizeType length{static_cast<SizeType>(batchSize * nbatches)};
    VecType y{length, sunctx};
    Vec2D y2d((y.View()).data(), nbatches, batchSize);

    Kokkos::parallel_for(
      "fill_y", Kokkos::RangePolicy<ExecSpace>(0, nbatches),
      KOKKOS_LAMBDA(const SizeType i) {
        y2d(i, 0) = u0;
        y2d(i, 1) = v0;
        y2d(i, 2) = w0;
      });

    // Create vector of absolute tolerances
    VecType abstol{length, sunctx};
    N_VConst(SUN_RCONST(1.0e-10), abstol);

    // Create CVODE using Backward Differentiation Formula methods
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }

    // Initialize the integrator and set the ODE right-hand side function
    int retval = CVodeInit(cvode_mem, f, T0, y);
    if (check_flag(retval, "CVodeInit")) { return 1; }

    // Attach the user data structure
    retval = CVodeSetUserData(cvode_mem, &udata);
    if (check_flag(retval, "CVodeSetUserData")) { return 1; }

    // Specify the scalar relative tolerance and vector absolute tolerances
    retval = CVodeSVtolerances(cvode_mem, SUN_RCONST(1.0e-6), abstol);
    if (check_flag(retval, "CVodeSVtolerances")) { return 1; }

    // Create the matrix and linear solver objects
    std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;
    std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

    if (solver_type == 0)
    {
      // Create Kokkos dense block diagonal matrix
      A = std::make_unique<MatType>(nbatches, batchSize, batchSize, sunctx);

      // Create Kokkos batched dense linear solver
      LS = std::make_unique<LSType>(sunctx);

      // Attach the matrix and linear solver to CVODE
      retval = CVodeSetLinearSolver(cvode_mem, LS->get(), A->get());
      if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }

      // Set the user-supplied Jacobian function
      retval = CVodeSetJacFn(cvode_mem, Jac);
      if (check_flag(retval, "CVodeSetJacFn")) { return 1; }
    }
    else
    {
      // Create matrix-free GMRES linear solver
      LS = std::make_unique<sundials::experimental::SUNLinearSolverView>(
        SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx));

      // Attach the linear solver to CVODE
      retval = CVodeSetLinearSolver(cvode_mem, LS->get(), nullptr);
      if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }
    }

    // Final time and time between outputs
    const sunrealtype Tf    = SUN_RCONST(10.0);
    const sunrealtype dTout = SUN_RCONST(1.0);

    // Number of output times
    const int Nt = static_cast<int>(ceil(Tf / dTout));

    // Current time and first output time
    sunrealtype t    = T0;
    sunrealtype tout = T0 + dTout;

    // Initial output
    Vec2DHost y2d_h((y.HostView()).data(), nbatches, batchSize);
    sundials::kokkos::CopyFromDevice(y);
    Kokkos::fence();
    std::cout << "At t = " << t << std::endl;
    for (int j = 0; j < nbatches; j += 10)
    {
      std::cout << "  batch " << j << ": y = " << y2d_h(j, 0) << " "
                << y2d_h(j, 1) << " " << y2d_h(j, 2) << std::endl;
    }

    // Loop over output times
    for (int iout = 0; iout < Nt; iout++)
    {
      // Advance in time
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      if (check_flag(retval, "CVode")) { break; }

      // Output solution from some batches
      sundials::kokkos::CopyFromDevice(y);
      Kokkos::fence();
      std::cout << "At t = " << t << std::endl;
      for (int j = 0; j < nbatches; j += 10)
      {
        std::cout << "  batch " << j << ": y = " << y2d_h(j, 0) << " "
                  << y2d_h(j, 1) << " " << y2d_h(j, 2) << std::endl;
      }

      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

    // Print some final statistics
    long int nst, nfe, nsetups, nje, nni, ncfn, netf;

    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_flag(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_flag(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_flag(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_flag(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_flag(retval, "CVodeGetNumNonlinSolvConvFails");
    retval = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_flag(retval, "CVodeGetNumJacEvals");

    std::cout << "\nFinal Statistics:\n"
              << "  Steps            = " << nst << "\n"
              << "  RHS evals        = " << nfe << "\n"
              << "  LS setups        = " << nsetups << "\n"
              << "  Jac evals        = " << nje << "\n"
              << "  NLS iters        = " << nni << "\n"
              << "  NLS fails        = " << ncfn << "\n"
              << "  Error test fails = " << netf << "\n";

    // Free objects
    CVodeFree(&cvode_mem);
  }
  Kokkos::finalize();

  return 0;
}

/* -----------------------------------------------------------------------------
 * User-supplied functions called by CVODE
 * ---------------------------------------------------------------------------*/

// Right hand side function dy/dt = f(t,y)
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto nbatches  = udata->nbatches;
  const auto batchSize = udata->batchSize;

  const auto a  = udata->a;
  const auto b  = udata->b;
  const auto ep = udata->ep;

  Vec2D y2d(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  Vec2D ydot2d(N_VGetDeviceArrayPointer(ydot), nbatches, batchSize);

  Kokkos::parallel_for(
    "RHS", Kokkos::RangePolicy<ExecSpace>(0, nbatches),
    KOKKOS_LAMBDA(const SizeType i) {
      auto u       = y2d(i, 0);
      auto v       = y2d(i, 1);
      auto w       = y2d(i, 2);
      ydot2d(i, 0) = a - (w + ONE) * u + v * u * u;
      ydot2d(i, 1) = w * u - v * u * u;
      ydot2d(i, 2) = (b - w) / ep - w * u;
    });

  return 0;
}

// Jacobian of f(t,y)
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto udata  = static_cast<UserData*>(user_data);
  auto y_data = sundials::kokkos::GetVec<VecType>(y)->View();
  auto J_data = sundials::kokkos::GetDenseMat<MatType>(J)->View();

  const auto nbatches  = udata->nbatches;
  const auto batchSize = udata->batchSize;

  const auto ep = udata->ep;
  Vec2D y2d(N_VGetDeviceArrayPointer(y), nbatches, batchSize);

  Kokkos::parallel_for(
    "Jac", Kokkos::RangePolicy<ExecSpace>(0, nbatches),
    KOKKOS_LAMBDA(const SizeType i) {
      // get y values
      auto u = y2d(i, 0);
      auto v = y2d(i, 1);
      auto w = y2d(i, 2);

      // first col of block
      J_data(i, 0, 0) = -(w + ONE) + TWO * u * v;
      J_data(i, 1, 0) = u * u;
      J_data(i, 2, 0) = -u;

      // second col of block
      J_data(i, 0, 1) = u * u;
      J_data(i, 1, 1) = -u * u;
      J_data(i, 2, 1) = u;

      // third col of block
      J_data(i, 0, 2) = -w;
      J_data(i, 1, 2) = ZERO;
      J_data(i, 2, 2) = -ONE / ep - u;
    });

  return 0;
}
