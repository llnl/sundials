/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * ARKODE main for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_lsrkstep.h"
#include "diffusion_2D.hpp"
#include "sunadaptcontroller/sunadaptcontroller_imexgus.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

struct UserOptions
{
  // Integrator settings
  sunrealtype rtol       = SUN_RCONST(1.0e-5);  // relative tolerance
  sunrealtype atol       = SUN_RCONST(1.0e-10); // absolute tolerance
  sunrealtype hfixed     = ZERO;                // fixed step size
  int order              = 3;                   // ARKode method order
  std::string controller = "I";                 // step size adaptivity method
  int maxsteps           = 0;                   // max steps between outputs
  int onestep            = 0;    // one step mode, number of steps
  bool linear            = true; // linearly implicit RHS
  bool implicit = true; // implicit (ARKStep) vs explicit STS (LSRKStep)
  ARKODE_LSRKMethodType lsrkmethod = ARKODE_LSRK_RKC_2; // LSRK method type

  // Linear solver and preconditioner settings
  std::string ls       = "cg";  // linear solver to use
  bool preconditioning = true;  // preconditioner on/off
  bool lsinfo          = false; // output residual history
  int liniters         = 20;    // number of linear iterations
  int msbp             = 0;     // preconditioner setup frequency
  sunrealtype epslin   = ZERO;  // linear solver tolerance factor

  // Helper functions
  int parse_args(vector<string>& args, bool outproc);
  void help();
  void print();
};

// -----------------------------------------------------------------------------
// LSRKStep-specific dominant eigenvalue function prototype
// -----------------------------------------------------------------------------

static int dom_eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
                   sunrealtype* lambdaI, void* user_data, N_Vector temp1,
                   N_Vector temp2, N_Vector temp3);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // Reusable error-checking flag
  int flag;

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) { return 1; }

  // Create SUNDIALS context
  MPI_Comm comm    = MPI_COMM_WORLD;
  SUNContext ctx   = nullptr;
  SUNProfiler prof = nullptr;

  flag = SUNContext_Create(comm, &ctx);
  if (check_flag(&flag, "SUNContextCreate", 1)) { return 1; }

  flag = SUNContext_GetProfiler(ctx, &prof);
  if (check_flag(&flag, "SUNContext_GetProfiler", 1)) { return 1; }

  // Add scope so objects are destroyed before MPI_Finalize
  {
    SUNDIALS_CXX_MARK_FUNCTION(prof);

    // MPI process ID
    int myid;
    flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (check_flag(&flag, "MPI_Comm_rank", 1)) { return 1; }

    bool outproc = (myid == 0);

    // --------------------------
    // Parse command line inputs
    // --------------------------

    UserData udata(prof);
    UserOptions uopts;
    UserOutput uout;

    vector<string> args(argv + 1, argv + argc);

    flag = udata.parse_args(args, outproc);
    if (check_flag(&flag, "UserData::parse_args", 1)) { return 1; }

    flag = uopts.parse_args(args, outproc);
    if (check_flag(&flag, "UserOptions::parse_args", 1)) { return 1; }

    flag = uout.parse_args(args, outproc);
    if (check_flag(&flag, "UserOutput::parse_args", 1)) { return 1; }

    // Check for unparsed inputs
    if (args.size() > 0)
    {
      if (find(args.begin(), args.end(), "--help") == args.end())
      {
        cerr << "ERROR: Unknown inputs: ";
        for (auto i = args.begin(); i != args.end(); ++i) { cerr << *i << ' '; }
        cerr << endl;
      }
      return 1;
    }

    // Return with error on unsupported LSRK method type
    if (!uopts.implicit)
    {
      if ((uopts.lsrkmethod != ARKODE_LSRK_RKC_2) &&
          (uopts.lsrkmethod != ARKODE_LSRK_RKL_2))
      {
        cerr << "ERROR: illegal lsrkmethod" << endl;
        return 1;
      }
    }

    // -----------------------------
    // Setup parallel decomposition
    // -----------------------------

    flag = udata.setup();
    if (check_flag(&flag, "UserData::setup", 1)) { return 1; }

    if (outproc)
    {
      udata.print();
      uopts.print();
      uout.print();
    }

    // ---------------
    // Create vectors
    // ---------------

    // Create vector for solution
#if defined(USE_HIP)
    N_Vector u = N_VMake_MPIPlusX(udata.comm_c,
                                  N_VNew_Hip(udata.nodes_loc, ctx), ctx);
    if (check_flag((void*)u, "N_VMake_MPIPlusX", 0)) return 1;
#elif defined(USE_CUDA)
    N_Vector u = N_VMake_MPIPlusX(udata.comm_c,
                                  N_VNew_Cuda(udata.nodes_loc, ctx), ctx);
    if (check_flag((void*)u, "N_VMake_MPIPlusX", 0)) return 1;
#else
    N_Vector u = N_VNew_Parallel(udata.comm_c, udata.nodes_loc, udata.nodes, ctx);
    if (check_flag((void*)u, "N_VNew_Parallel", 0)) { return 1; }
#endif

    // Set initial condition
    flag = Solution(ZERO, u, &udata);
    if (check_flag(&flag, "Solution", 1)) { return 1; }

    // Create vector for error
    if (udata.forcing)
    {
      uout.error = N_VClone(u);
      if (check_flag((void*)(uout.error), "N_VClone", 0)) { return 1; }
    }

    // Set up implicit solver, if applicable
    SUNLinearSolver LS = nullptr;
    SUNMatrix A        = nullptr;
#if defined(USE_SUPERLU_DIST)
    // SuperLU-DIST objects
    SuperMatrix A_super;
    gridinfo_t grid;
    dLUstruct_t A_lu;
    dScalePermstruct_t A_scaleperm;
    dSOLVEstruct_t A_solve;
    SuperLUStat_t A_stat;
    superlu_dist_options_t A_opts;
    sunrealtype* A_data      = nullptr;
    sunindextype* A_col_idxs = nullptr;
    sunindextype* A_row_ptrs = nullptr;
#endif
    if (uopts.implicit)
    {
      // ---------------------
      // Create linear solver
      // ---------------------

      // Create linear solver

      int prectype = (uopts.preconditioning) ? SUN_PREC_RIGHT : SUN_PREC_NONE;

      if (uopts.ls == "cg")
      {
        LS = SUNLinSol_PCG(u, prectype, uopts.liniters, ctx);
        if (check_flag((void*)LS, "SUNLinSol_PCG", 0)) { return 1; }
      }
      else if (uopts.ls == "gmres")
      {
        LS = SUNLinSol_SPGMR(u, prectype, uopts.liniters, ctx);
        if (check_flag((void*)LS, "SUNLinSol_SPGMR", 0)) { return 1; }
      }
      else
      {
#if defined(USE_SUPERLU_DIST)
        // Initialize SuperLU-DIST grid
        superlu_gridinit(udata.comm_c, udata.npx, udata.npy, &grid);

        // Create arrays for CSR matrix: data, column indices, and row pointers
        sunindextype nnz_loc = 5 * udata.nodes_loc;

        A_data = (sunrealtype*)malloc(nnz_loc * sizeof(sunrealtype));
        if (check_flag((void*)A_data, "malloc Adata", 0)) return 1;

        A_col_idxs = (sunindextype*)malloc(nnz_loc * sizeof(sunindextype));
        if (check_flag((void*)A_col_idxs, "malloc Acolind", 0)) return 1;

        A_row_ptrs =
          (sunindextype*)malloc((udata.nodes_loc + 1) * sizeof(sunindextype));
        if (check_flag((void*)A_row_ptrs, "malloc Arowptr", 0)) return 1;

        // Create and initialize SuperLU_DIST structures
        dCreate_CompRowLoc_Matrix_dist(&A_super, udata.nodes, udata.nodes,
                                       nnz_loc, udata.nodes_loc, 0, A_data,
                                       A_col_idxs, A_row_ptrs, SLU_NR_loc,
                                       SLU_D, SLU_GE);
        dScalePermstructInit(udata.nodes, udata.nodes, &A_scaleperm);
        dLUstructInit(udata.nodes, &A_lu);
        PStatInit(&A_stat);
        set_default_options_dist(&A_opts);
        A_opts.PrintStat = NO;

        // SUNDIALS structures
        A = SUNMatrix_SLUNRloc(&A_super, &grid, ctx);
        if (check_flag((void*)A, "SUNMatrix_SLUNRloc", 0)) return 1;

        LS = SUNLinSol_SuperLUDIST(u, A, &grid, &A_lu, &A_scaleperm, &A_solve,
                                   &A_stat, &A_opts, ctx);
        if (check_flag((void*)LS, "SUNLinSol_SuperLUDIST", 0)) return 1;

        uopts.preconditioning = false;
#else
        std::cerr
          << "ERROR: Benchmark was not built with SuperLU_DIST enabled\n";
        return 1;
#endif
      }

      // Allocate preconditioner workspace
      if (uopts.preconditioning)
      {
        udata.diag = N_VClone(u);
        if (check_flag((void*)(udata.diag), "N_VClone", 0)) { return 1; }
      }
    }

    // ----------------------
    // Setup ARKStep/LSRKStep
    // ----------------------

    // Create integrator
    void* arkode_mem = nullptr;
    if (uopts.implicit)
    {
      arkode_mem = ARKStepCreate(nullptr, diffusion, ZERO, u, ctx);
      if (check_flag((void*)arkode_mem, "ARKStepCreate", 0)) { return 1; }
    }
    else
    {
      arkode_mem = LSRKStepCreateSTS(diffusion, ZERO, u, ctx);
      if (check_flag((void*)arkode_mem, "LSRKStepCreateSTS", 0)) { return 1; }
    }

    // Specify tolerances
    flag = ARKodeSStolerances(arkode_mem, uopts.rtol, uopts.atol);
    if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

    // Attach user data
    flag = ARKodeSetUserData(arkode_mem, (void*)&udata);
    if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

    // Configure implicit solver
    if (uopts.implicit)
    {
      // Attach linear solver
      flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
      if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }

#if defined(USE_SUPERLU_DIST)
      if (uopts.ls == "sludist")
      {
        ARKodeSetJacFn(arkode_mem, diffusion_jac);
        if (check_flag(&flag, "ARKodeSetJacFn", 1)) return 1;
      }
#endif

      if (uopts.preconditioning)
      {
        // Attach preconditioner
        flag = ARKodeSetPreconditioner(arkode_mem, PSetup, PSolve);
        if (check_flag(&flag, "ARKodeSetPreconditioner", 1)) { return 1; }

        // Set linear solver setup frequency (update preconditioner)
        flag = ARKodeSetLSetupFrequency(arkode_mem, uopts.msbp);
        if (check_flag(&flag, "ARKodeSetLSetupFrequency", 1)) { return 1; }
      }

      // Set linear solver tolerance factor
      flag = ARKodeSetEpsLin(arkode_mem, uopts.epslin);
      if (check_flag(&flag, "ARKodeSetEpsLin", 1)) { return 1; }

      // Select method order
      flag = ARKodeSetOrder(arkode_mem, uopts.order);
      if (check_flag(&flag, "ARKodeSetOrder", 1)) { return 1; }

      // Specify linearly implicit non-time-dependent RHS
      if (uopts.linear)
      {
        flag = ARKodeSetLinear(arkode_mem, 0);
        if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }
      }
    }
    else // Configure explicit STS solver
    {
      // Select LSRK method
      flag = LSRKStepSetSTSMethod(arkode_mem, uopts.lsrkmethod);
      if (check_flag(&flag, "LSRKStepSetSTSMethod", 1)) { return 1; }

      // Provide dominant eigenvalue function
      flag = LSRKStepSetDomEigFn(arkode_mem, dom_eig);
      if (check_flag(&flag, "LSRKStepSetDomEigFn", 1)) { return 1; }
    }

    // Set fixed step size or adaptivity method
    if (uopts.hfixed > ZERO)
    {
      flag = ARKodeSetFixedStep(arkode_mem, uopts.hfixed);
      if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
    }
    else
    {
      flag = ARKodeSetAdaptControllerByName(arkode_mem, uopts.controller.c_str());
      if (check_flag(&flag, "ARKodeSetAdaptControllerByName", 1)) { return 1; }
    }

    // Set max steps between outputs
    flag = ARKodeSetMaxNumSteps(arkode_mem, uopts.maxsteps);
    if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

    // Set stopping time
    flag = ARKodeSetStopTime(arkode_mem, udata.tf);
    if (check_flag(&flag, "ARKodeSetStopTime", 1)) { return 1; }

    // -----------------------
    // Loop over output times
    // -----------------------

    // Optionally run in one step mode for a fixed number of time steps (helpful
    // for debugging)
    int stepmode = ARK_NORMAL;

    if (uopts.onestep)
    {
      uout.nout = uopts.onestep;
      stepmode  = ARK_ONE_STEP;
    }

    sunrealtype t     = ZERO;
    sunrealtype dTout = udata.tf / uout.nout;
    sunrealtype tout  = dTout;

    // Initial output
    flag = uout.open(&udata);
    if (check_flag(&flag, "UserOutput::open", 1)) { return 1; }

    flag = uout.write(t, u, &udata);
    if (check_flag(&flag, "UserOutput::write", 1)) { return 1; }

    for (int iout = 0; iout < uout.nout; iout++)
    {
      SUNDIALS_MARK_BEGIN(prof, "Evolve");

      // Evolve in time
      flag = ARKodeEvolve(arkode_mem, tout, u, &t, stepmode);
      if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }

      SUNDIALS_MARK_END(prof, "Evolve");

      // Output solution and error
      flag = uout.write(t, u, &udata);
      if (check_flag(&flag, "UserOutput::write", 1)) { return 1; }

      // Update output time
      tout += dTout;
      tout = (tout > udata.tf) ? udata.tf : tout;
    }

    // Close output
    flag = uout.close(&udata);
    if (check_flag(&flag, "UserOutput::close", 1)) { return 1; }

    // --------------
    // Final outputs
    // --------------

    // Print final integrator stats
    if (outproc)
    {
      cout << "Final integrator statistics:" << endl;
      flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
      if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }
    }

    // ---------
    // Clean up
    // ---------

    // Free MPI Cartesian communicator
    MPI_Comm_free(&(udata.comm_c));

    ARKodeFree(&arkode_mem);
    if (uopts.implicit)
    {
      SUNLinSolFree(LS);

      // Free the SuperLU_DIST structures (also frees user allocated arrays
      // A_data, A_col_idxs, and A_row_ptrs)
#if defined(USE_SUPERLU_DIST)
      if (uopts.ls == "sludist")
      {
        PStatFree(&A_stat);
        dScalePermstructFree(&A_scaleperm);
        dLUstructFree(&A_lu);
        Destroy_CompRowLoc_Matrix_dist(&A_super);
        superlu_gridexit(&grid);
      }
#endif
    }

    // Free vectors
#if defined(USE_HIP) || defined(USE_CUDA)
    N_VDestroy(N_VGetLocalVector_MPIPlusX(u));
#endif
    N_VDestroy(u);
  }
  // Close scope so objects are destroyed before MPI_Finalize

  SUNContext_Free(&ctx);

  // Finalize MPI
  flag = MPI_Finalize();
  return 0;
}

// -----------------------------------------------------------------------------
// Dominant eigenvalue estimation function
// -----------------------------------------------------------------------------

static int dom_eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
                   sunrealtype* lambdaI, void* user_data, N_Vector temp1,
                   N_Vector temp2, N_Vector temp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Fill in spectral radius value
  *lambdaR = -SUN_RCONST(8.0) * std::max(udata->kx / udata->dx / udata->dx,
                                         udata->ky / udata->dy / udata->dy);
  *lambdaI = SUN_RCONST(0.0);

  // return with success
  return 0;
}

// -----------------------------------------------------------------------------
// UserOptions Helper functions
// -----------------------------------------------------------------------------

int UserOptions::parse_args(vector<string>& args, bool outproc)
{
  vector<string>::iterator it;

  it = find(args.begin(), args.end(), "--help");
  if (it != args.end())
  {
    if (outproc) { help(); }
    return 0;
  }

  it = find(args.begin(), args.end(), "--rtol");
  if (it != args.end())
  {
    rtol = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--atol");
  if (it != args.end())
  {
    atol = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--fixedstep");
  if (it != args.end())
  {
    hfixed = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--explicitSTS");
  if (it != args.end())
  {
    implicit = false;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--order");
  if (it != args.end())
  {
    order = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--lsrkmethod");
  if (it != args.end())
  {
    lsrkmethod = (ARKODE_LSRKMethodType)stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--controller");
  if (it != args.end())
  {
    controller = *(it + 1);
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--maxsteps");
  if (it != args.end())
  {
    maxsteps = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--onestep");
  if (it != args.end())
  {
    onestep = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--nonlinear");
  if (it != args.end())
  {
    linear = false;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--ls");
  if (it != args.end())
  {
    ls = *(it + 1);
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--lsinfo");
  if (it != args.end())
  {
    lsinfo = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--liniters");
  if (it != args.end())
  {
    liniters = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--msbp");
  if (it != args.end())
  {
    msbp = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--epslin");
  if (it != args.end())
  {
    epslin = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  return 0;
}

// Print command line options
void UserOptions::help()
{
  cout << endl;
  cout << "Integrator command line options:" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --atol <atol>           : absolute tolerance" << endl;
  cout << "  --controller <ctr>      : time step adaptivity controller" << endl;
  cout << "  --fixedstep <step>      : used fixed step size" << endl;
  cout << "  --explicitSTS           : use LSRKStep (instead of ARKStep)" << endl;
  cout << endl;
  cout << "Implicit (ARKStep) solver command line options:" << endl;
  cout << "  --nonlinear             : disable linearly implicit flag" << endl;
  cout << "  --order <ord>           : method order" << endl;
  cout << "  --ls <cg|gmres|sludist> : linear solver" << endl;
  cout << "  --lsinfo                : output residual history" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --noprec                : disable preconditioner" << endl;
  cout << "  --msbp <steps>          : max steps between prec setups" << endl;
  cout << endl;
  cout << "Explicit STS (LSRKStep) solver command line options:" << endl;
  cout << "  --lsrkmethod            : LSRK method choice" << endl;
}

// Print user options
void UserOptions::print()
{
  cout << endl;
  cout << " Integrator options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " rtol        = " << rtol << endl;
  cout << " atol        = " << atol << endl;
  cout << " controller  = " << controller << endl;
  cout << " hfixed      = " << hfixed << endl;
  cout << " --------------------------------- " << endl;
  cout << endl;
  if (implicit)
  {
    cout << " ARKStep options:" << endl;
    cout << " --------------------------------- " << endl;
    cout << " order       = " << order << endl;
    cout << " max steps   = " << maxsteps << endl;
    cout << " linear RHS  = " << linear << endl;
    cout << " --------------------------------- " << endl;
    cout << endl;
    if (ls == "sludist")
    {
      cout << " Linear solver options:" << endl;
      cout << " --------------------------------- " << endl;
#if defined(HAVE_HIP)
      cout << " LS       = SuperLU_DIST (HIP enabled)" << endl;
#elif defined(HAVE_CUDA)
      cout << " LS       = SuperLU_DIST (CUDA enabled)" << endl;
#else
      cout << " LS       = SuperLU_DIST" << endl;
#endif
      cout << " LS info  = " << lsinfo << endl;
      cout << " msbp     = " << msbp << endl;
      cout << " --------------------------------- " << endl;
    }
    else
    {
      cout << " Linear solver options:" << endl;
      cout << " --------------------------------- " << endl;
      cout << " LS       = " << ls << endl;
      cout << " precond  = " << preconditioning << endl;
      cout << " LS info  = " << lsinfo << endl;
      cout << " LS iters = " << liniters << endl;
      cout << " msbp     = " << msbp << endl;
      cout << " epslin   = " << epslin << endl;
      cout << " --------------------------------- " << endl;
    }
  }
  else
  {
    cout << " LSRKStep options:" << endl;
    cout << " --------------------------------- " << endl;
    switch (lsrkmethod)
    {
    case (ARKODE_LSRK_RKC_2): cout << " method = RKC_2 " << endl; break;
    case (ARKODE_LSRK_RKL_2): cout << " method = RKL_2 " << endl; break;
    default: cout << " ERROR: illegal lsrkmethod " << endl;
    }
    cout << " --------------------------------- " << endl;
  }
}

//---- end of file ----
