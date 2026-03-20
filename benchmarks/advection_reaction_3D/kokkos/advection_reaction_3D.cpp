/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                David J. Gardner, Cody J. Balos @ LLNL
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
 * This benchmark problem simulates the advection and reaction of three
 * chemical species, u, v, and w, in a three dimensional domain. The reaction
 * mechanism is a variation of the Brusselator problem from chemical kinetics.
 * This is a PDE system with 3 components, Y = [u,v,w], satisfying the
 * equations,
 *
 *    u_t = -c * dot(grad,u) + A - (w+1) * u + v * u^2
 *    v_t = -c * dot(grad,v) + w * u - v * u^2
 *    w_t = -c * dot(grad,w) + (B - w) / ep - w * u
 *
 * for t in [0,tf], X = (x,y,z) where in (x,y,z) in [0,xmax] with periodic
 * boundary conditions. The initial condition is
 *
 *    u(0,X) = k1 * A / k4 + p(X)
 *    v(0,X) = k2 * k4 * B / (k1 * k3 * A) + p(X)
 *    w(0,X) = 3.0 + p(X)
 *    p(X)   = alpha * e^( -((X - mu)^T Sigma^{-1} (x-mu)) / (2*sqrt(|Sigma|*(2pi)^3)) )
 *
 * alpha = 0.1, mu = (xmax/2.0, xmax/2.0, xmax/2.0), and Sigma = diag(xmax/4.0).
 * The reaction rates are set so k_1 = k_2 = k_3 = k_4 = k, and k_5 = k_6
 * = 1/5e-6. The spatial derivatives are discretized with first-order upwind
 * finite differences. NOUT outputs are printed at equal intervals, and run
 * statistics are printed at the end.
 *
 * Command line options:
 *   --help             prints this message
 *   --dont-save        do not save the solution to the filesystem at the nout interval (default is to save)
 *   --output-dir       the directory where all output files will be written
 *   --nout <int>       number of output times
 *   --method           ERK, ARK-DIRK, ARK-IMEX (default), CV-BDF, CV-ADAMS, IDA
 *   --nls              nonlinear solver to use; options are newton,
 *                      tl-newton (task-local newton), or fixedpoint
 *   --fpaccel          the number of fixed-point acceleration vectors to use
 *                      (only valid when using fixedpoint nonlinear solver)
 *   --nopre            turn off preconditioning
 *   --order <int>      the method order to use
 *   --npts <int>       number of mesh points in each direction
 *   --xmax <sunrealtype>  maximum value of x (size of domain)
 *   --tf <sunrealtype>    final time
 *   --A <sunrealtype>     A parameter value
 *   --B <sunrealtype>     B parameter value
 *   --k <sunrealtype>     reaction rate
 *   --c <sunrealtype>     advection speed
 *   --rtol <sunrealtype>  relative tolerance
 *   --atol <sunrealtype>  absolute tolerance
 * --------------------------------------------------------------------------*/

#include "advection_reaction_3D.hpp"

/* Main Program */
int main(int argc, char* argv[])
{
  SUNContext ctx;

  /* Initialize MPI */
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  /* Create SUNDIALS context */
  SUNContext_Create(comm, &ctx);

  /* Initialize Kokkos */
  Kokkos::initialize(argc, argv);
  {
    /* General problem variables */
    N_Vector y = NULL;   /* empty solution vector        */
    UserData udata(ctx); /* user data                    */
    UserOptions uopt;    /* user options                 */
    int retval;          /* reusable error-checking flag */

    SUNDIALS_CXX_MARK_FUNCTION(udata.prof);

    /* Process input arguments and set up the problem */
    retval = SetupProblem(argc, argv, &udata, &uopt, ctx);
    if (check_retval(&retval, "SetupProblem", 1, udata.myid))
    {
      MPI_Abort(comm, 1);
    }

    /* Create solution vector (on-node and MPI-parallel versions) */
    SUNVector yloc{(unsigned int)udata.grid->neq, ctx};
    y = N_VMake_MPIPlusX(udata.comm, yloc, ctx);
    if (check_retval((void*)y, "N_VMake_MPIPlusX", 0, udata.myid))
    {
      MPI_Abort(comm, 1);
    }

    /* Set the initial condition */
    retval = SetIC(y, &udata);
    if (check_retval(&retval, "SetIC", 1, udata.myid)) { MPI_Abort(comm, 1); }

    /* Output spatial mesh to disk (add extra point for periodic BC) */
    if (udata.myid == 0 && uopt.nout > 0)
    {
      char fname[MXSTR];
      snprintf(fname, MXSTR, "%s/mesh.txt", uopt.outputdir.c_str());
      udata.grid->MeshToFile(fname);
    }

    /* Integrate in time */
    if (uopt.method == "ERK")
    {
      retval = EvolveProblemExplicit(y, &udata, &uopt);
    }
    else if (uopt.method == "ARK-DIRK")
    {
      retval = EvolveProblemDIRK(y, &udata, &uopt);
    }
    else if (uopt.method == "ARK-IMEX")
    {
      retval = EvolveProblemIMEX(y, &udata, &uopt);
    }
    else if (uopt.method == "CV-BDF")
    {
      retval = EvolveProblemBDF(y, &udata, &uopt);
    }
    else if (uopt.method == "CV-ADAMS")
    {
      retval = EvolveProblemAdams(y, &udata, &uopt);
    }
    else if (uopt.method == "IDA")
    {
      retval = EvolveDAEProblem(y, &udata, &uopt);
    }
    if (check_retval(&retval, "Evolve", 1, udata.myid)) { MPI_Abort(comm, 1); }

    /* Clean up */
    N_VDestroy(y);
  }
  Kokkos::finalize();
  SUNContext_Free(&ctx);
  MPI_Finalize();
  return (0);
}

/* Destructor for problem data */
UserData::~UserData()
{
  /* close output streams */
  if (uopt->nout > 0)
  {
    if (UFID) { fclose(UFID); }
    if (VFID) { fclose(VFID); }
    if (WFID) { fclose(WFID); }
    if (TFID && myid == 0) { fclose(TFID); }
  }

  /* free solution masks */
  if (umask != nullptr)
  {
    N_VDestroy(umask);
    umask = nullptr;
  }
  if (vmask != nullptr)
  {
    N_VDestroy(vmask);
    vmask = nullptr;
  }
  if (wmask != nullptr)
  {
    N_VDestroy(wmask);
    wmask = nullptr;
  }

  /* free the parallel grid */
  delete grid;
}

/* --------------------------------------------------------------
 * Communication functions
 * --------------------------------------------------------------*/

/* Fills send buffers before exchanging neighbor information */
int FillSendBuffers(N_Vector y, UserData* udata)
{
  /* Shortcuts */
  const sunrealtype c = udata->c;
  const int nxl       = udata->grid->nxl;
  const int nyl       = udata->grid->nyl;
  const int nzl       = udata->grid->nzl;
  const int dof       = udata->grid->dof;

  /* Create 4D view of the vector */
  Vec4D Yview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y)), nxl, nyl,
              nzl, dof);

  if (c > 0.0)
  {
    /* Flow moving in the positive directions uses backward difference. */

    /* Create 4D views of send buffers */
    Vec4D Esend(udata->grid->GetSendView("EAST"), 1, nyl, nzl, dof);
    Vec4D Nsend(udata->grid->GetSendView("NORTH"), nxl, 1, nzl, dof);
    Vec4D Fsend(udata->grid->GetSendView("FRONT"), nxl, nyl, 1, dof);

    /* Fill buffers on device */
    Kokkos::parallel_for(
      "FillEastBuffer", Range3D({0, 0, 0}, {nyl, nzl, dof}),
      KOKKOS_LAMBDA(int j, int k, int l) {
        Esend(0, j, k, l) = Yview(nxl - 1, j, k, l);
      });
    Kokkos::parallel_for(
      "FillNorthBuffer", Range3D({0, 0, 0}, {nxl, nzl, dof}),
      KOKKOS_LAMBDA(int i, int k, int l) {
        Nsend(i, 0, k, l) = Yview(i, nyl - 1, k, l);
      });
    Kokkos::parallel_for(
      "FillFrontBuffer", Range3D({0, 0, 0}, {nxl, nyl, dof}),
      KOKKOS_LAMBDA(int i, int j, int l) {
        Fsend(i, j, 0, l) = Yview(i, j, nzl - 1, l);
      });
  }
  else if (c < 0.0)
  {
    /* Flow moving in the negative directions uses forward difference. */

    /* Create 4D views of send buffers */
    Vec4D Wsend(udata->grid->GetSendView("WEST"), 1, nyl, nzl, dof);
    Vec4D Ssend(udata->grid->GetSendView("SOUTH"), nxl, 1, nzl, dof);
    Vec4D Bsend(udata->grid->GetSendView("BACK"), nxl, nyl, 1, dof);

    /* Fill buffers on device */
    Kokkos::parallel_for(
      "FillWestBuffer", Range3D({0, 0, 0}, {nyl, nzl, dof}),
      KOKKOS_LAMBDA(int j, int k, int l) {
        Wsend(0, j, k, l) = Yview(0, j, k, l);
      });
    Kokkos::parallel_for(
      "FillSouthBuffer", Range3D({0, 0, 0}, {nxl, nzl, dof}),
      KOKKOS_LAMBDA(int i, int k, int l) {
        Ssend(i, 0, k, l) = Yview(i, 0, k, l);
      });
    Kokkos::parallel_for(
      "FillBackBuffer", Range3D({0, 0, 0}, {nxl, nyl, dof}),
      KOKKOS_LAMBDA(int i, int j, int l) {
        Bsend(i, j, 0, l) = Yview(i, j, 0, l);
      });
  }

  return (0);
}

/* --------------------------------------------------------------
 * Problem setup
 * --------------------------------------------------------------*/

/* Parses the CLI arguments */
static int ParseArgs(int argc, char* argv[], UserData* udata, UserOptions* uopt)
{
  /* check for input args */
  if (argc > 1)
  {
    /* loop over input args and get value */
    for (int i = 1; i < argc; i++)
    {
      string argvi(argv[i]);

      if (argvi.compare("--help") == 0)
      {
        InputError(argv[0]);
        return (-1);
      }
      else if (argvi.compare("--nout") == 0) { uopt->nout = atoi(argv[++i]); }
      else if (argvi.compare("--dont-save") == 0) { uopt->save = 0; }
      else if (argvi.compare("--output-dir") == 0)
      {
        if (strlen(argv[i + 1]) > MXSTR)
        {
          if (udata->myid == 0)
          {
            fprintf(stderr, "ERROR: output directory string is too long\n");
          }
          return (-1);
        }
        uopt->outputdir = argv[++i];
      }
      else if (argvi.compare("--npts") == 0) { uopt->npts = atoi(argv[++i]); }
      else if (argvi.compare("--npxyz") == 0)
      {
        uopt->npxyz[0] = atoi(argv[++i]);
        uopt->npxyz[1] = atoi(argv[++i]);
        uopt->npxyz[2] = atoi(argv[++i]);
      }
      else if (argvi.compare("--xmax") == 0)
      {
        udata->xmax = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--A") == 0)
      {
        udata->A = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--B") == 0)
      {
        udata->B = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--k") == 0)
      {
        udata->k1 = strtod(argv[++i], NULL);
        udata->k2 = strtod(argv[++i], NULL);
        udata->k3 = strtod(argv[++i], NULL);
        udata->k4 = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--c") == 0)
      {
        udata->c = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--order") == 0) { uopt->order = atoi(argv[++i]); }
      else if (argvi.compare("--method") == 0)
      {
        uopt->method = string(argv[++i]);
        if (uopt->method != "ERK" && uopt->method != "ARK-DIRK" &&
            uopt->method != "ARK-IMEX" && uopt->method != "CV-BDF" &&
            uopt->method != "CV-ADAMS" && uopt->method != "IDA")
        {
          fprintf(stderr, "ERROR: unknown method\n");
          InputError(argv[0]);
          return (-1);
        }
      }
      else if (argvi.compare("--fpaccel") == 0)
      {
        uopt->fpaccel = atoi(argv[++i]);
      }
      else if (argvi.compare("--nls") == 0)
      {
        uopt->nls = string(argv[++i]);
        if (uopt->nls != "newton" && uopt->nls != "tl-newton" &&
            uopt->nls != "fixedpoint" && uopt->nls != "none")
        {
          fprintf(stderr, "ERROR: unknown nls\n");
          InputError(argv[0]);
          return (-1);
        }
      }
      else if (argvi.compare("--nopre") == 0) { uopt->precond = 0; }
      else if (argvi.compare("--fused") == 0) { uopt->fused = 1; }
      else if (argvi.compare("--tf") == 0)
      {
        uopt->tf = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--rtol") == 0)
      {
        uopt->rtol = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--atol") == 0)
      {
        uopt->atol = strtod(argv[++i], NULL);
      }
      else
      {
        InputError(argv[0]);
        return (-1);
      }
    }
  }

  /* Explicit method uses no nonlinear solver */
  if (uopt->method == "ERK") { uopt->nls = "none"; }

  /* CV Adams method only uses fixedpoint nonlinear solver */
  if (uopt->method == "CV-ADAMS") { uopt->nls = "fixedpoint"; }

  return (0);
}

/* Fills the mask vector for the component so that
   u = y .* umask, v = y .* vmask, w = y .* wmask */
int ComponentMask(N_Vector mask, const int component, const UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* Shortcuts */
  const int nxl = udata->grid->nxl;
  const int nyl = udata->grid->nyl;
  const int nzl = udata->grid->nzl;
  const int dof = udata->grid->dof;

  /* Create 4D view of mask data */
  Vec4D maskview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(mask)),
                 nxl, nyl, nzl, dof);

  /* Fill mask data */
  N_VConst(0.0, mask);
  Kokkos::parallel_for(
    "Fill_mask", Range3D({0, 0, 0}, {nxl, nyl, nzl}),
    KOKKOS_LAMBDA(int i, int j, int k) { maskview(i, j, k, component) = 1.0; });

  return 0;
}

/* Parses the CLI arguments and sets up the problem */
int SetupProblem(int argc, char* argv[], UserData* udata, UserOptions* uopt,
                 SUNContext ctx)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* MPI variables */
  udata->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(udata->comm, &udata->myid);
  MPI_Comm_size(udata->comm, &udata->nprocs);

  /* Default problem parameters */
  udata->add_reactions = true;
  udata->xmax          = 1.0;
  udata->A             = 1.0;
  udata->B             = 3.5;
  udata->k1            = 1.0;
  udata->k2            = 1.0;
  udata->k3            = 1.0;
  udata->k4            = 1.0;
  udata->k5            = 1.0 / 5.0e-6;
  udata->k6            = 1.0 / 5.0e-6;
  udata->c             = 0.01;
  udata->uopt          = uopt;
  udata->TFID          = NULL;
  udata->UFID          = NULL;
  udata->VFID          = NULL;
  udata->WFID          = NULL;
  udata->nnlfi         = 0;

  /* Set default integrator options */
  uopt->npxyz[0]  = 0;          /* number of processesors in x */
  uopt->npxyz[1]  = 0;          /* number of processesors in y */
  uopt->npxyz[2]  = 0;          /* number of processesors in z */
  uopt->npts      = 100;        /* number of mesh points in each direction */
  uopt->order     = 3;          /* method order             */
  uopt->method    = "ARK-DIRK"; /* stepper/method           */
  uopt->t0        = 0.0;        /* initial time             */
  uopt->tf        = 10.0;       /* final time               */
  uopt->rtol      = 1.0e-6;     /* relative tolerance       */
  uopt->atol      = 1.0e-9;     /* absolute tolerance       */
  uopt->nls       = "newton";   /* default to newton, when appropriate */
  uopt->fpaccel   = 3;  /* default number of fixed point acceleration vectors */
  uopt->precond   = 1;  /* by default, precondition when appropriate */
  uopt->fused     = 0;  /* use fused vector ops     */
  uopt->save      = 1;  /* save solution to disk    */
  uopt->nout      = 10; /* number of output times   */
  uopt->outputdir = "."; /* output directory         */

  /* Parse CLI args and set udata/uopt appropriately */
  int retval = ParseArgs(argc, argv, udata, uopt);
  if (check_retval((void*)&retval, "ParseArgs", 1, udata->myid)) { return -1; }

  /* Setup the parallel decomposition */
  const sunindextype npts[] = {uopt->npts, uopt->npts, uopt->npts};
  const sunrealtype amax[]  = {0.0, 0.0, 0.0};
  const sunrealtype bmax[]  = {udata->xmax, udata->xmax, udata->xmax};
  udata->grid = new ParallelGrid<sunindextype>(&udata->comm, amax, bmax, npts,
                                               3, BoundaryType::PERIODIC,
                                               StencilType::UPWIND, udata->c,
                                               uopt->npxyz);

  /* Create the solution masks */
  SUNVector* umaskloc = new SUNVector((unsigned int)udata->grid->neq, ctx);
  udata->umask        = N_VMake_MPIPlusX(udata->comm, *umaskloc, ctx);
  if (check_retval((void*)udata->umask, "N_VMake_MPIPlusX", 0, udata->myid))
  {
    MPI_Abort(udata->comm, 1);
  }
  SUNVector* vmaskloc = new SUNVector((unsigned int)udata->grid->neq, ctx);
  udata->vmask        = N_VMake_MPIPlusX(udata->comm, *vmaskloc, ctx);
  if (check_retval((void*)udata->vmask, "N_VMake_MPIPlusX", 0, udata->myid))
  {
    MPI_Abort(udata->comm, 1);
  }
  SUNVector* wmaskloc = new SUNVector((unsigned int)udata->grid->neq, ctx);
  udata->wmask        = N_VMake_MPIPlusX(udata->comm, *wmaskloc, ctx);
  if (check_retval((void*)udata->wmask, "N_VMake_MPIPlusX", 0, udata->myid))
  {
    MPI_Abort(udata->comm, 1);
  }
  ComponentMask(udata->umask, 0, udata);
  ComponentMask(udata->vmask, 1, udata);
  ComponentMask(udata->wmask, 2, udata);

  /* Open output files for results */
  if (uopt->save)
  {
    char fname[MXSTR];
    if (udata->myid == 0)
    {
      sprintf(fname, "%s/t.%06d.txt", uopt->outputdir.c_str(), udata->myid);
      udata->TFID = fopen(fname, "w");
    }

    sprintf(fname, "%s/u.%06d.txt", uopt->outputdir.c_str(), udata->myid);
    udata->UFID = fopen(fname, "w");

    sprintf(fname, "%s/v.%06d.txt", uopt->outputdir.c_str(), udata->myid);
    udata->VFID = fopen(fname, "w");

    sprintf(fname, "%s/w.%06d.txt", uopt->outputdir.c_str(), udata->myid);
    udata->WFID = fopen(fname, "w");
  }

  /* Print problem setup */
  if (udata->myid == 0)
  {
    printf("\n\t\tAdvection-Reaction Test Problem\n\n");
    printf("Using the MPI+Kokkos NVECTOR");
#if defined(USE_CUDA)
    printf(" with the CUDA back-end\n");
#elif defined(USE_HIP)
    printf(" with the HIP back-end\n");
#elif defined(USE_OPENMP)
    printf(" with the OpenMP back-end and %i threads\n", omp_get_max_threads());
#else
    printf(" with the serial back-end\n");
#endif
    printf("Number of Processors = %li\n", (long int)udata->nprocs);
    udata->grid->PrintInfo();
    printf("Problem Parameters:\n");
    printf("  A = %g\n", udata->A);
    printf("  B = %g\n", udata->B);
    printf("  k = %g\n", udata->k1);
    printf("  c = %g\n", udata->c);
    printf("Integrator Options:\n");
    printf("  order            = %d\n", uopt->order);
    printf("  method           = %s\n", uopt->method.c_str());
    printf("  nonlinear solver = %s\n", uopt->nls.c_str());
    printf("  fpaccel          = %d\n", uopt->fpaccel);
    printf("  preconditioner   = %d\n", uopt->precond);
    printf("  fused vector ops = %d\n", uopt->fused);
    printf("  t0               = %g\n", uopt->t0);
    printf("  tf               = %g\n", uopt->tf);
    printf("  reltol           = %.1e\n", uopt->rtol);
    printf("  abstol           = %.1e\n", uopt->atol);
    printf("  nout             = %d\n", uopt->nout);
    printf("Output directory: %s\n", uopt->outputdir.c_str());
  }

  /* return success */
  return (0);
}

/* Compute the 3D Gaussian function. */
KOKKOS_FUNCTION
static void Gaussian3D(sunrealtype& x, sunrealtype& y, sunrealtype& z,
                       sunrealtype xmax)
{
  /* Gaussian distribution defaults */
  const sunrealtype alpha   = 0.1;
  const sunrealtype mu[]    = {xmax / SUN_RCONST(2.0), xmax / SUN_RCONST(2.0),
                               xmax / SUN_RCONST(2.0)};
  const sunrealtype sigma[] = {xmax / SUN_RCONST(4.0), xmax / SUN_RCONST(4.0),
                               xmax / SUN_RCONST(4.0)}; // Sigma = diag(sigma)

  /* denominator = 2*sqrt(|Sigma|*(2pi)^3) */
  const sunrealtype denom =
    2.0 * sqrt((sigma[0] * sigma[1] * sigma[2]) * pow(2 * M_PI, 3));
  x = alpha * exp(-((x - mu[0]) * (x - mu[0]) * (1.0 / sigma[0])) / denom);
  y = alpha * exp(-((y - mu[1]) * (y - mu[1]) * (1.0 / sigma[1])) / denom);
  z = alpha * exp(-((z - mu[2]) * (z - mu[2]) * (1.0 / sigma[2])) / denom);
}

/* Initial condition function */
int SetIC(N_Vector yvec, UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* Variable shortcuts */
  const int nxl          = udata->grid->nxl;
  const int nyl          = udata->grid->nyl;
  const int nzl          = udata->grid->nzl;
  const int dof          = udata->grid->dof;
  const sunrealtype dx   = udata->grid->dx;
  const sunrealtype dy   = udata->grid->dy;
  const sunrealtype dz   = udata->grid->dz;
  const sunrealtype xmax = udata->xmax;
  const sunrealtype A    = udata->A;
  const sunrealtype B    = udata->B;
  const sunrealtype k1   = udata->k1;
  const sunrealtype k2   = udata->k2;
  const sunrealtype k3   = udata->k3;
  const sunrealtype k4   = udata->k4;
  const int xcrd         = udata->grid->coords[0];
  const int ycrd         = udata->grid->coords[1];
  const int zcrd         = udata->grid->coords[2];

  /* Steady state solution */
  const sunrealtype us = k1 * A / k4;
  const sunrealtype vs = k2 * k4 * B / (k1 * k3 * A);
  const sunrealtype ws = 3.0;

  /* Create 4D view of y */
  Vec4D yview(N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(yvec)), nxl,
              nyl, nzl, dof);

  /* Gaussian perturbation of the steady state solution */
  Kokkos::parallel_for(
    "SetIC", Range3D({0, 0, 0}, {nxl, nyl, nzl}),
    KOKKOS_LAMBDA(int i, int j, int k) {
      sunrealtype x = (xcrd * nxl + i) * dx;
      sunrealtype y = (ycrd * nyl + j) * dy;
      sunrealtype z = (zcrd * nzl + k) * dz;
      Gaussian3D(x, y, z, xmax);
      const sunrealtype p = x + y + z;
      yview(i, j, k, 0)   = us + p;
      yview(i, j, k, 1)   = vs + p;
      yview(i, j, k, 2)   = ws + p;
    });

  /* Return success */
  return (0);
}

/* Write time and solution to disk */
int WriteOutput(sunrealtype t, N_Vector y, UserData* udata, UserOptions* uopt)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* output current solution norm to screen */
  sunrealtype N = (sunrealtype)udata->grid->npts();
  sunrealtype u = N_VWL2Norm(y, udata->umask);
  u             = sqrt(u * u / N);
  sunrealtype v = N_VWL2Norm(y, udata->vmask);
  v             = sqrt(v * v / N);
  sunrealtype w = N_VWL2Norm(y, udata->wmask);
  w             = sqrt(w * w / N);
  if (udata->myid == 0)
  {
    printf("     %10.6f   %10.6f   %10.6f   %10.6f\n", t, u, v, w);
    std::fflush(stdout);
  }

  if (uopt->save)
  {
    /* Copy solution data to host mirror view */
    SUNVector* ylocal =
      sundials::kokkos::GetVec<SUNVector>(N_VGetLocalVector_MPIPlusX(y));
    sundials::kokkos::CopyFromDevice(*ylocal);

    /* output the times to disk */
    if (udata->myid == 0 && udata->TFID)
    {
      fprintf(udata->TFID, " %.16e\n", t);
      std::fflush(udata->TFID);
    }

    /* create 4D view of host data */
    const int nxl = udata->grid->nxl;
    const int nyl = udata->grid->nyl;
    const int nzl = udata->grid->nzl;
    const int dof = udata->grid->dof;
    Vec4DHost yview(N_VGetArrayPointer(N_VGetLocalVector_MPIPlusX(y)), nxl, nyl,
                    nzl, dof);

    /* output results to disk */
    for (int i = 0; i < nxl; i++)
    {
      for (int j = 0; j < nyl; j++)
      {
        for (int k = 0; k < nzl; k++)
        {
          fprintf(udata->UFID, " %.16e", yview(i, j, k, 0));
          fprintf(udata->VFID, " %.16e", yview(i, j, k, 1));
          fprintf(udata->WFID, " %.16e", yview(i, j, k, 2));
        }
      }
    }

    fprintf(udata->UFID, "\n");
    fprintf(udata->VFID, "\n");
    fprintf(udata->WFID, "\n");
    std::fflush(udata->UFID);
    std::fflush(udata->VFID);
    std::fflush(udata->WFID);
  }

  return (0);
}

void InputError(char* name)
{
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0)
  {
    fprintf(stderr, "\nERROR: Invalid command line input\n");
    fprintf(stderr, "\nCommand line options for %s\n", name);
    fprintf(stderr, "  --help                    prints this message\n");
    fprintf(stderr, "  --output-dir              the directory where all "
                    "output files will be written (default is the CWD)\n");
    fprintf(stderr, "  --nout <int>              number of output times to "
                    "print (default is 10)\n");
    fprintf(stderr,
            "  --dont-save               do not save the solution to the "
            "filesystem at the nout interval (default is to save)\n");
    fprintf(stderr, "  --method                  ERK, ARK-DIRK, ARK-IMEX "
                    "(default), CV-BDF, CV-ADAMS, IDA\n");
    fprintf(stderr, "  --fpaccel                 the number of fixed-point "
                    "acceleration vectors to use (only valid when using "
                    "fixedpoint nonlinear solver)\n");
    fprintf(stderr, "  --nls                     nonlinear solver to use "
                    "(newton, tl-newton (task-local newton), fixedpoint)\n");
    fprintf(stderr, "  --nopre                   do not precondition the "
                    "linear system\n");
    fprintf(stderr, "  --order <int>             the method order to use\n");
    fprintf(stderr, "  --npts <int>              number of mesh points in each "
                    "direction\n");
    fprintf(stderr, "  --npxyz <int> <int> <int> number of processors in each "
                    "direction (0 forces MPI to decide)\n");
    fprintf(stderr, "  --xmax <sunrealtype>         maximum value of x (size "
                    "of domain)\n");
    fprintf(stderr, "  --tf <sunrealtype>           final time\n");
    fprintf(stderr, "  --A <sunrealtype>            A parameter value\n");
    fprintf(stderr, "  --B <sunrealtype>            B parameter value\n");
    fprintf(stderr, "  --k <sunrealtype>            reaction rate\n");
    fprintf(stderr, "  --c <sunrealtype>            advection speed\n");
    fprintf(stderr, "  --rtol <sunrealtype>         relative tolerance\n");
    fprintf(stderr, "  --atol <sunrealtype>         absolute tolerance\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
