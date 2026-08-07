/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * ---------------------------------------------------------------------------*/

#ifndef ADVECTION_REACTION_3D_HPP
#define ADVECTION_REACTION_3D_HPP

#include <RAJA/RAJA.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <nvector/nvector_mpiplusx.h>
#include <sundials/sundials_core.hpp>

#include "ParallelGrid.hpp"
#include "backends.hpp"
#include "check_retval.h"

using std::string;
using sundials_tools::BoundaryType;
using sundials_tools::ParallelGrid;
using sundials_tools::StencilType;

/* Maximum size of output directory string */
constexpr int MXSTR = 2048;

/*
 * Data structure for problem options
 */

struct UserOptions
{
  int npxyz[3];      /* number of processors in x,y,z */
  sunindextype npts; /* number of spatial mesh points */
  sunrealtype t0;    /* initial time                  */
  sunrealtype tf;    /* final time                    */
  sunrealtype rtol;  /* relative tolerance            */
  sunrealtype atol;  /* absolute tolerance            */
  int order;         /* method order                  */
  string method;     /* method string                 */
  string nls;        /* nonlinear solver to use       */
  int fpaccel;       /* number of fixedpoint vectors  */
  int precond;       /* to precondition or not        */
  int fused;         /* use fused vector ops          */
  int nout;          /* number of outputs             */
  int save;          /* save solution to disk         */
  const char* outputdir;
};

/*
 * Data structure for problem specific data
 */

struct UserData
{
  SUNContext ctx;
  SUNProfiler prof;

  /* MPI data */
  MPI_Comm comm;
  int myid;
  int nprocs;
  MPI_Request req[2];

  /* should reactions be added to the advection or not */
  bool add_reactions;

  /* file handles for output */
  FILE* TFID; /* time output file pointer     */
  FILE* UFID; /* solution output file pointer */
  FILE* VFID;
  FILE* WFID;

  /* solution masks */
  N_Vector umask;
  N_Vector vmask;
  N_Vector wmask;

  /* problem parameters */
  sunrealtype xmax; /* maximum x value              */
  sunrealtype A;    /* concentration of species A   */
  sunrealtype B;    /* w source rate                */
  sunrealtype k1;   /* reaction rates               */
  sunrealtype k2;
  sunrealtype k3;
  sunrealtype k4;
  sunrealtype k5;
  sunrealtype k6;
  sunrealtype c; /* advection coefficient        */

  /* parallel mesh */
  ParallelGrid<sunrealtype, sunindextype>* grid;

  /* count of implicit function evals by the task local nonlinear solver */
  long int nnlfi;

  /* integrator options */
  UserOptions* uopt;

  /* constructor that takes the context */
  UserData(SUNContext ctx)
    : ctx(ctx),
      umask(nullptr),
      vmask(nullptr),
      wmask(nullptr),
      uopt(nullptr),
      TFID(nullptr),
      UFID(nullptr),
      VFID(nullptr),
      WFID(nullptr)
  {
    SUNContext_GetProfiler(ctx, &prof);
  }

  /* destructor frees the problem data */
  ~UserData();
};

/*
 * Functions to evolve the solution (defined by the drivers)
 */

/* function that does ARKStep setup and evolves the solution with a DIRK method */
extern int EvolveProblemDIRK(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does ARKStep setup and evolves the solution with an IMEX method */
extern int EvolveProblemIMEX(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does ERKStep setup and evolves the solution */
extern int EvolveProblemExplicit(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does CVODE BDF setup and evolves the solution */
extern int EvolveProblemBDF(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does CVODE Adams setup and evolves the solution */
extern int EvolveProblemAdams(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does IDA BDF setup and evolves the solution */
extern int EvolveDAEProblem(N_Vector y, UserData* udata, UserOptions* uopt);

/*
 * Helper functions
 */

/* function to set initial condition */
int SetIC(N_Vector y, UserData* udata);

/* function to fill neighbor data */
int FillSendBuffers(N_Vector y, UserData* udata);

/* functions for processing command line args */
int SetupProblem(int argc, char* argv[], UserData* udata, UserOptions* uopt,
                 SUNMemoryHelper memhelper, SUNContext ctx);
void InputError(char* name);
int ComponentMask(N_Vector mask, const int component, const UserData* udata);

/* function to write solution to disk */
int WriteOutput(sunrealtype t, N_Vector y, UserData* udata, UserOptions* uopt);

#endif
