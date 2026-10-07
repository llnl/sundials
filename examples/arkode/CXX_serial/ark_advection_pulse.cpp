/* -----------------------------------------------------------------------------
 * Programmer(s): SUNDIALS development team
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
 * Example problem:
 *
 * This example solves the one-dimensional wave/advection equation
 *
 *   u_t = u_x,  x in [0, 1],
 *
 * with the initial condition u(0, x) = 0. Since the wave speed is -1, data
 * enters through the right boundary. A smooth time-dependent pulse is imposed
 * at x = 1 and advects left through the domain.
 *
 * The method-of-lines ODE is posed as the DAE residual u' - f(t,u) = 0 and
 * solved with IDA or with PDAEStep using IDA on each partition.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <example_utilities.hpp>

#include <arkode/arkode.h>
#include <arkode/arkode_pdaestep.h>
#include <ida/ida.h>
#include <ida/ida_ls.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

struct UserData
{
  sunindextype n;
  sunrealtype dx;
  sunrealtype pulse_start;
  sunrealtype pulse_width;
};

struct Options
{
  std::string solver;
  int partitions;
  int ida_max_order;
  sunrealtype pdae_h;
};

struct PdaeUserData
{
  UserData model;
  int partitions;
  sunindextype coupling_size;
  std::vector<sunindextype> starts;
  std::vector<sunindextype> sizes;
};

struct SolutionStats
{
  sunrealtype max_norm;
  sunrealtype l1_norm;
  sunrealtype x_peak;
};

static constexpr sunrealtype ZERO = SUN_RCONST(0.0);
static constexpr sunrealtype ONE  = SUN_RCONST(1.0);

static sunrealtype Pulse(sunrealtype t, const UserData& data);
static sunrealtype AdvectionRhsValue(sunrealtype uright, sunrealtype u,
                                     const UserData& data);
static void AdvectionRhs(sunrealtype t, N_Vector u, N_Vector udot,
                         const UserData& data);
static void PartitionLayout(sunindextype n, int partitions,
                            std::vector<sunindextype>& starts,
                            std::vector<sunindextype>& sizes);
static int PdaeComponentResidualForPartition(int partition, sunrealtype t,
                                             N_Vector y, N_Vector w, N_Vector yp,
                                             N_Vector res, void* user_data);
static int PdaeComponentJacTimesForPartition(int partition, sunrealtype t,
                                             N_Vector y, N_Vector w,
                                             N_Vector yp, N_Vector res,
                                             N_Vector v, N_Vector Jv,
                                             sunrealtype cj, void* user_data,
                                             N_Vector tmp1, N_Vector tmp2);
static int PdaeAlgebraicResidual(sunrealtype t, N_Vector y, N_Vector w,
                                 N_Vector res, void* user_data);
static int PdaeAlgebraicJacTimes(sunrealtype t, N_Vector y, N_Vector w,
                                 N_Vector v, N_Vector Jv, void* user_data,
                                 N_Vector tmp);
static void CopyPdaeSolution(N_Vector y, const PdaeUserData& data, N_Vector u);
static int Residual(sunrealtype t, N_Vector u, N_Vector up, N_Vector res,
                    void* user_data);
static int JacTimes(sunrealtype t, N_Vector u, N_Vector up, N_Vector res,
                    N_Vector v, N_Vector Jv, sunrealtype cj, void* user_data,
                    N_Vector tmp1, N_Vector tmp2);
static SolutionStats GetStats(N_Vector u, const UserData& data);
static void PrintHeader();
static void PrintSolution(int step, sunrealtype t, N_Vector u,
                          const UserData& data);
static int WriteSolution(std::ofstream& outfile, sunrealtype t, N_Vector u,
                         const UserData& data);
static int PrintL1Error(N_Vector u, N_Vector uref);
static int PrintIdaFinalStats(void* ida_mem);
static int PrintPdaeFinalStats(void* arkode_mem,
                               const std::vector<void*>& ida_mems);
static int RunIda(const UserData& data, sunrealtype t0, sunrealtype tf,
                  sunrealtype reltol, sunrealtype abstol, int nout,
                  SUNContext ctx, N_Vector ufinal = NULL,
                  bool print_output = true, int max_order = 1);
static int RunPdae(const UserData& data, const Options& opts, sunrealtype t0,
                   sunrealtype tf, sunrealtype reltol, sunrealtype abstol,
                   int nout, N_Vector uref, SUNContext ctx);
static int ParseArgs(std::vector<std::string>& args, UserData& data,
                     sunrealtype& tf, sunrealtype& reltol, sunrealtype& abstol,
                     int& nout, Options& opts);

int main(int argc, char* argv[])
{
  sundials::Context ctx;

  UserData data;
  data.n           = 2001;
  data.dx          = ONE / static_cast<sunrealtype>(data.n);
  data.pulse_start = SUN_RCONST(0.05);
  data.pulse_width = SUN_RCONST(0.20);

  sunrealtype t0     = ZERO;
  sunrealtype tf     = SUN_RCONST(1.20);
  sunrealtype reltol = SUN_RCONST(1.0e-4);
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  int nout           = 12;

  Options opts;
  opts.solver        = "ida";
  opts.partitions    = 5;
  opts.ida_max_order = 1;
  opts.pdae_h        = SUN_RCONST(1.0e-4);

  std::vector<std::string> args(argv + 1, argv + argc);
  int flag = ParseArgs(args, data, tf, reltol, abstol, nout, opts);
  if (flag != 0) { return flag; }

  N_Vector uref = N_VNew_Serial(data.n, ctx);
  if (check_ptr(uref, "N_VNew_Serial")) { return 1; }

  const sunrealtype ref_reltol = SUN_RCONST(1.0e-8);
  const sunrealtype ref_abstol = SUN_RCONST(1.0e-12);
  flag = RunIda(data, t0, tf, ref_reltol, ref_abstol, nout, ctx, uref, false, 0);

  if (flag == 0 && opts.solver == "pdae")
  {
    flag = RunPdae(data, opts, t0, tf, reltol, abstol, nout, uref, ctx);
  }
  else if (flag == 0)
  {
    N_Vector u = N_VNew_Serial(data.n, ctx);
    if (check_ptr(u, "N_VNew_Serial")) { flag = 1; }
    else
    {
      flag = RunIda(data, t0, tf, reltol, abstol, nout, ctx, u, true,
                    opts.ida_max_order);
      if (flag == 0) { flag = PrintL1Error(u, uref); }
      N_VDestroy(u);
    }
  }

  N_VDestroy(uref);

  return flag;
}

static int RunIda(const UserData& data, sunrealtype t0, sunrealtype tf,
                  sunrealtype reltol, sunrealtype abstol, int nout,
                  SUNContext ctx, N_Vector ufinal, bool print_output,
                  int max_order)
{
  if (print_output)
  {
    std::cout << "\nAdvection pulse problem solved with IDA:\n";
    std::cout << "   PDE          = u_t = u_x\n";
    std::cout << "   mesh points  = " << data.n << "\n";
    std::cout << "   dx           = " << data.dx << "\n";
    std::cout << "   t0           = " << t0 << "\n";
    std::cout << "   tf           = " << tf << "\n";
    std::cout << "   pulse start  = " << data.pulse_start << "\n";
    std::cout << "   pulse width  = " << data.pulse_width << "\n";
    std::cout << "   reltol       = " << reltol << "\n";
    std::cout << "   abstol       = " << abstol << "\n\n";
  }

  N_Vector u = N_VNew_Serial(data.n, ctx);
  if (check_ptr(u, "N_VNew_Serial")) { return 1; }

  N_Vector up = N_VClone(u);
  if (check_ptr(up, "N_VClone")) { return 1; }

  N_VConst(ZERO, u);
  AdvectionRhs(t0, u, up, data);

  void* ida_mem = IDACreate(ctx);
  if (check_ptr(ida_mem, "IDACreate")) { return 1; }

  int flag = IDAInit(ida_mem, Residual, t0, u, up);
  if (check_flag(flag, "IDAInit")) { return 1; }

  flag = IDASetUserData(ida_mem, const_cast<UserData*>(&data));
  if (check_flag(flag, "IDASetUserData")) { return 1; }

  if (max_order > 0)
  {
    flag = IDASetMaxOrd(ida_mem, max_order);
    if (check_flag(flag, "IDASetMaxOrd")) { return 1; }
  }

  flag = IDASStolerances(ida_mem, reltol, abstol);
  if (check_flag(flag, "IDASStolerances")) { return 1; }

  flag = IDASetMaxNumSteps(ida_mem, -1);
  if (check_flag(flag, "IDASetMaxNumSteps")) { return 1; }

  SUNLinearSolver LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 10, ctx);
  if (check_ptr(LS, "SUNLinSol_SPGMR")) { return 1; }

  flag = IDASetLinearSolver(ida_mem, LS, NULL);
  if (check_flag(flag, "IDASetLinearSolver")) { return 1; }

  flag = IDASetJacTimes(ida_mem, NULL, JacTimes);
  if (check_flag(flag, "IDASetJacTimes")) { return 1; }

  std::ofstream outfile;
  if (print_output)
  {
    outfile.open("ark_advection_pulse.txt");
    outfile << "# vars: t x u\n";
    outfile << std::scientific;
    outfile << std::setprecision(std::numeric_limits<sunrealtype>::digits10);

    PrintHeader();
    PrintSolution(0, t0, u, data);
    flag = WriteSolution(outfile, t0, u, data);
    if (check_flag(flag, "WriteSolution")) { return 1; }
  }

  sunrealtype tret = t0;
  for (int iout = 1; iout <= nout; iout++)
  {
    sunrealtype tout = tf * static_cast<sunrealtype>(iout) /
                       static_cast<sunrealtype>(nout);

    flag = IDASetStopTime(ida_mem, tout);
    if (check_flag(flag, "IDASetStopTime")) { return 1; }

    flag = IDASolve(ida_mem, tout, &tret, u, up, IDA_NORMAL);
    if (check_flag(flag, "IDASolve")) { return 1; }

    if (print_output)
    {
      PrintSolution(iout, tret, u, data);
      flag = WriteSolution(outfile, tret, u, data);
      if (check_flag(flag, "WriteSolution")) { return 1; }
    }
  }

  if (ufinal) { N_VScale(ONE, u, ufinal); }

  if (print_output)
  {
    std::cout << "\n";
    flag = PrintIdaFinalStats(ida_mem);
  }
  else
  {
    flag = 0;
  }

  if (outfile.is_open()) { outfile.close(); }

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  N_VDestroy(up);
  N_VDestroy(u);

  return flag;
}

static sunrealtype Pulse(sunrealtype t, const UserData& data)
{
  if (t < data.pulse_start || t > data.pulse_start + data.pulse_width)
  {
    return ZERO;
  }

  const sunrealtype pi = SUN_RCONST(3.141592653589793238462643383279502884);
  const sunrealtype s  = (t - data.pulse_start) / data.pulse_width;
  const sunrealtype p  = std::sin(pi * s);

  return p * p;
}

static sunrealtype AdvectionRhsValue(sunrealtype uright, sunrealtype u,
                                     const UserData& data)
{
  return (uright - u) / data.dx;
}

static void AdvectionRhs(sunrealtype t, N_Vector u, N_Vector udot,
                         const UserData& data)
{
  sunrealtype* udata    = N_VGetArrayPointer(u);
  sunrealtype* udotdata = N_VGetArrayPointer(udot);

  for (sunindextype i = 0; i < data.n; i++)
  {
    sunrealtype uright = (i == data.n - 1) ? Pulse(t, data) : udata[i + 1];
    udotdata[i]        = AdvectionRhsValue(uright, udata[i], data);
  }
}

static void PartitionLayout(sunindextype n, int partitions,
                            std::vector<sunindextype>& starts,
                            std::vector<sunindextype>& sizes)
{
  starts.resize(partitions);
  sizes.resize(partitions);

  sunindextype first = 0;
  sunindextype base  = n / partitions;
  sunindextype extra = n % partitions;

  for (int p = 0; p < partitions; p++)
  {
    sizes[p]  = base + ((p < extra) ? 1 : 0);
    starts[p] = first;
    first += sizes[p];
  }
}

static int PdaeComponentResidualForPartition(int partition, sunrealtype t,
                                             N_Vector y, N_Vector w, N_Vector yp,
                                             N_Vector res, void* user_data)
{
  PdaeUserData* data = static_cast<PdaeUserData*>(user_data);
  if (partition < 0 || partition >= data->partitions) { return -1; }

  N_Vector x_vec    = N_VGetSubvector_ManyVector(y, 0);
  N_Vector z_vec    = N_VGetSubvector_ManyVector(y, 1);
  N_Vector xp_vec   = N_VGetSubvector_ManyVector(yp, 0);
  N_Vector xres_vec = N_VGetSubvector_ManyVector(res, 0);
  N_Vector zres_vec = N_VGetSubvector_ManyVector(res, 1);

  sunrealtype* x    = N_VGetArrayPointer(x_vec);
  sunrealtype* z    = N_VGetArrayPointer(z_vec);
  sunrealtype* xp   = N_VGetArrayPointer(xp_vec);
  sunrealtype* xres = N_VGetArrayPointer(xres_vec);
  sunrealtype* zres = N_VGetArrayPointer(zres_vec);
  sunrealtype* wptr = N_VGetArrayPointer(w);

  const sunindextype nlocal = data->sizes[partition];
  for (sunindextype i = 0; i < nlocal; i++)
  {
    sunrealtype uright = (i == nlocal - 1) ? z[0] : x[i + 1];
    xres[i]            = xp[i] - AdvectionRhsValue(uright, x[i], data->model);
  }

  zres[0] = z[0] - ((partition == data->partitions - 1) ? Pulse(t, data->model)
                                                        : wptr[partition]);

  return 0;
}

static int PdaeAlgebraicResidual(sunrealtype t, N_Vector y, N_Vector w,
                                 N_Vector res, void* user_data)
{
  PdaeUserData* data = static_cast<PdaeUserData*>(user_data);
  sunrealtype* wptr  = N_VGetArrayPointer(w);

  for (int p = 0; p < data->partitions; p++)
  {
    N_Vector z_vec    = PDAEStepGetAlgebraicSubvector(y, p);
    N_Vector zres_vec = N_VGetSubvector_ManyVector(res, p);
    sunrealtype* z    = N_VGetArrayPointer(z_vec);
    sunrealtype* zres = N_VGetArrayPointer(zres_vec);

    zres[0] = z[0] -
              ((p == data->partitions - 1) ? Pulse(t, data->model) : wptr[p]);
  }

  N_Vector wres_vec = N_VGetSubvector_ManyVector(res, data->partitions);
  sunrealtype* wres = N_VGetArrayPointer(wres_vec);
  for (sunindextype i = 0; i < data->coupling_size; i++) { wres[i] = wptr[i]; }

  for (int p = 0; p < data->partitions - 1; p++)
  {
    N_Vector x_vec = PDAEStepGetDifferentialSubvector(y, p + 1);
    sunrealtype* x = N_VGetArrayPointer(x_vec);
    wres[p]        = wptr[p] - x[0];
  }

  return 0;
}

static int PdaeAlgebraicJacTimes(sunrealtype t, N_Vector y, N_Vector w,
                                 N_Vector v, N_Vector Jv, void* user_data,
                                 N_Vector tmp)
{
  PdaeUserData* data = static_cast<PdaeUserData*>(user_data);

  (void)t;
  (void)y;
  (void)w;
  (void)tmp;

  N_Vector vw_vec = N_VGetSubvector_ManyVector(v, data->partitions);
  N_Vector Jw_vec = N_VGetSubvector_ManyVector(Jv, data->partitions);
  sunrealtype* vw = N_VGetArrayPointer(vw_vec);
  sunrealtype* Jw = N_VGetArrayPointer(Jw_vec);

  for (int p = 0; p < data->partitions; p++)
  {
    N_Vector vz_vec = N_VGetSubvector_ManyVector(v, p);
    N_Vector Jz_vec = N_VGetSubvector_ManyVector(Jv, p);
    sunrealtype* vz = N_VGetArrayPointer(vz_vec);
    sunrealtype* Jz = N_VGetArrayPointer(Jz_vec);

    Jz[0] = vz[0] - ((p == data->partitions - 1) ? ZERO : vw[p]);
  }

  for (sunindextype i = 0; i < data->coupling_size; i++) { Jw[i] = vw[i]; }

  return 0;
}

static int PdaeComponentJacTimesForPartition(int partition, sunrealtype t,
                                             N_Vector y, N_Vector w,
                                             N_Vector yp, N_Vector res,
                                             N_Vector v, N_Vector Jv,
                                             sunrealtype cj, void* user_data,
                                             N_Vector tmp1, N_Vector tmp2)
{
  PdaeUserData* data = static_cast<PdaeUserData*>(user_data);
  if (partition < 0 || partition >= data->partitions) { return -1; }

  (void)t;
  (void)w;
  (void)yp;
  (void)res;
  (void)tmp1;
  (void)tmp2;

  N_Vector vx_vec = N_VGetSubvector_ManyVector(v, 0);
  N_Vector vz_vec = N_VGetSubvector_ManyVector(v, 1);
  N_Vector Jx_vec = N_VGetSubvector_ManyVector(Jv, 0);
  N_Vector Jz_vec = N_VGetSubvector_ManyVector(Jv, 1);

  sunrealtype* vx = N_VGetArrayPointer(vx_vec);
  sunrealtype* vz = N_VGetArrayPointer(vz_vec);
  sunrealtype* Jx = N_VGetArrayPointer(Jx_vec);
  sunrealtype* Jz = N_VGetArrayPointer(Jz_vec);

  const sunrealtype invdx       = ONE / data->model.dx;
  const sunrealtype diag        = cj + invdx;
  const sunindextype nlocal     = data->sizes[partition];
  const sunindextype right_cell = nlocal - 1;

  for (sunindextype i = 0; i < nlocal; i++)
  {
    Jx[i] = diag * vx[i] - invdx * ((i == right_cell) ? vz[0] : vx[i + 1]);
  }

  Jz[0] = vz[0];

  return 0;
}

static void CopyPdaeSolution(N_Vector y, const PdaeUserData& data, N_Vector u)
{
  sunrealtype* udata = N_VGetArrayPointer(u);

  for (int p = 0; p < data.partitions; p++)
  {
    N_Vector x_vec = PDAEStepGetDifferentialSubvector(y, p);
    sunrealtype* x = N_VGetArrayPointer(x_vec);

    for (sunindextype i = 0; i < data.sizes[p]; i++)
    {
      udata[data.starts[p] + i] = x[i];
    }
  }
}

static int Residual(sunrealtype t, N_Vector u, N_Vector up, N_Vector res,
                    void* user_data)
{
  UserData* data = static_cast<UserData*>(user_data);
  AdvectionRhs(t, u, res, *data);
  N_VLinearSum(ONE, up, -ONE, res, res);
  return 0;
}

static int JacTimes(sunrealtype t, N_Vector u, N_Vector up, N_Vector res,
                    N_Vector v, N_Vector Jv, sunrealtype cj, void* user_data,
                    N_Vector tmp1, N_Vector tmp2)
{
  UserData* data     = static_cast<UserData*>(user_data);
  sunrealtype* vdata = N_VGetArrayPointer(v);
  sunrealtype* Jdata = N_VGetArrayPointer(Jv);

  (void)t;
  (void)u;
  (void)up;
  (void)res;
  (void)tmp1;
  (void)tmp2;

  const sunrealtype invdx = ONE / data->dx;
  const sunrealtype diag  = cj + invdx;

  for (sunindextype i = 0; i < data->n; i++)
  {
    Jdata[i] = diag * vdata[i] -
               ((i == data->n - 1) ? ZERO : invdx * vdata[i + 1]);
  }

  return 0;
}

static SolutionStats GetStats(N_Vector u, const UserData& data)
{
  const sunrealtype* udata = N_VGetArrayPointer(u);

  SolutionStats stats;
  stats.max_norm = ZERO;
  stats.l1_norm  = ZERO;
  stats.x_peak   = SUN_RCONST(0.5) * data.dx;

  sunindextype imax = 0;
  for (sunindextype i = 0; i < data.n; i++)
  {
    const sunrealtype abs_ui = std::abs(udata[i]);
    stats.l1_norm += data.dx * abs_ui;
    if (abs_ui > stats.max_norm)
    {
      stats.max_norm = abs_ui;
      imax           = i;
    }
  }

  stats.x_peak = (static_cast<sunrealtype>(imax) + SUN_RCONST(0.5)) * data.dx;

  return stats;
}

static void PrintHeader()
{
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << std::setw(8) << "output" << std::setw(22) << "t" << std::setw(22)
            << "max |u|" << std::setw(22) << "L1" << std::setw(22) << "x peak"
            << std::setw(22) << "u(1,t)"
            << "\n";
  for (int i = 0; i < 8 + 5 * 22; i++) { std::cout << "-"; }
  std::cout << "\n";
}

static void PrintSolution(int step, sunrealtype t, N_Vector u,
                          const UserData& data)
{
  SolutionStats stats = GetStats(u, data);
  std::cout << std::setw(8) << step << std::setw(22) << t << std::setw(22)
            << stats.max_norm << std::setw(22) << stats.l1_norm << std::setw(22)
            << stats.x_peak << std::setw(22) << Pulse(t, data) << "\n";
}

static int WriteSolution(std::ofstream& outfile, sunrealtype t, N_Vector u,
                         const UserData& data)
{
  if (!outfile.is_open()) { return -1; }

  const sunrealtype* udata = N_VGetArrayPointer(u);

  for (sunindextype i = 0; i < data.n; i++)
  {
    const sunrealtype x = (static_cast<sunrealtype>(i) + SUN_RCONST(0.5)) *
                          data.dx;
    outfile << t << " " << x << " " << udata[i] << std::endl;
  }
  outfile << std::endl;

  return 0;
}

static int PrintL1Error(N_Vector u, N_Vector uref)
{
  N_Vector err = N_VClone(u);
  if (check_ptr(err, "N_VClone")) { return 1; }

  N_VLinearSum(ONE, u, -ONE, uref, err);
  std::cout << "   Final solution 1-norm error = " << N_VL1Norm(err) << "\n";

  N_VDestroy(err);

  return 0;
}

static int PrintIdaFinalStats(void* ida_mem)
{
  long int nst     = 0;
  long int nre     = 0;
  long int nsetups = 0;
  long int nje     = 0;
  long int njv     = 0;
  long int nli     = 0;
  long int netf    = 0;
  long int nni     = 0;
  long int ncfn    = 0;

  int flag = IDAGetNumSteps(ida_mem, &nst);
  if (check_flag(flag, "IDAGetNumSteps")) { return 1; }

  flag = IDAGetNumResEvals(ida_mem, &nre);
  if (check_flag(flag, "IDAGetNumResEvals")) { return 1; }

  flag = IDAGetNumLinSolvSetups(ida_mem, &nsetups);
  if (check_flag(flag, "IDAGetNumLinSolvSetups")) { return 1; }

  flag = IDAGetNumJacEvals(ida_mem, &nje);
  if (check_flag(flag, "IDAGetNumJacEvals")) { return 1; }

  flag = IDAGetNumJtimesEvals(ida_mem, &njv);
  if (check_flag(flag, "IDAGetNumJtimesEvals")) { return 1; }

  flag = IDAGetNumLinIters(ida_mem, &nli);
  if (check_flag(flag, "IDAGetNumLinIters")) { return 1; }

  flag = IDAGetNumErrTestFails(ida_mem, &netf);
  if (check_flag(flag, "IDAGetNumErrTestFails")) { return 1; }

  flag = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  if (check_flag(flag, "IDAGetNumNonlinSolvIters")) { return 1; }

  flag = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  if (check_flag(flag, "IDAGetNumNonlinSolvConvFails")) { return 1; }

  std::cout << "Final Solver Statistics:\n";
  std::cout << "   Internal solver steps = " << nst << "\n";
  std::cout << "   Total residual evaluations = " << nre << "\n";
  std::cout << "   Total linear solver setups = " << nsetups << "\n";
  std::cout << "   Total Jacobian evaluations = " << nje << "\n";
  std::cout << "   Total J-times evaluations = " << njv << "\n";
  std::cout << "   Total linear iterations = " << nli << "\n";
  std::cout << "   Total error test failures = " << netf << "\n";
  std::cout << "   Total nonlinear iterations = " << nni << "\n";
  std::cout << "   Total nonlinear convergence failures = " << ncfn << "\n";

  return 0;
}

static int PrintPdaeFinalStats(void* arkode_mem,
                               const std::vector<void*>& ida_mems)
{
  long int nst          = 0;
  long int nlinsetups   = 0;
  long int nniters      = 0;
  long int nnfails      = 0;
  long int ida_nst      = 0;
  long int ida_nre      = 0;
  long int ida_nsetups  = 0;
  long int ida_njv      = 0;
  long int ida_nli      = 0;
  long int ida_netf     = 0;
  long int ida_nni      = 0;
  long int ida_ncfn     = 0;
  long int ida_nst_tmp  = 0;
  long int ida_nre_tmp  = 0;
  long int ida_nls_tmp  = 0;
  long int ida_njv_tmp  = 0;
  long int ida_nli_tmp  = 0;
  long int ida_netf_tmp = 0;
  long int ida_nni_tmp  = 0;
  long int ida_ncfn_tmp = 0;

  int flag = ARKodeGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "ARKodeGetNumSteps")) { return 1; }

  flag = PDAEStepGetNumLinSolvSetups(arkode_mem, &nlinsetups);
  if (check_flag(flag, "PDAEStepGetNumLinSolvSetups")) { return 1; }

  flag = PDAEStepGetNonlinSolvStats(arkode_mem, &nniters, &nnfails);
  if (check_flag(flag, "PDAEStepGetNonlinSolvStats")) { return 1; }

  for (void* ida_mem : ida_mems)
  {
    flag = IDAGetNumSteps(ida_mem, &ida_nst_tmp);
    if (check_flag(flag, "IDAGetNumSteps")) { return 1; }
    flag = IDAGetNumResEvals(ida_mem, &ida_nre_tmp);
    if (check_flag(flag, "IDAGetNumResEvals")) { return 1; }
    flag = IDAGetNumLinSolvSetups(ida_mem, &ida_nls_tmp);
    if (check_flag(flag, "IDAGetNumLinSolvSetups")) { return 1; }
    flag = IDAGetNumJtimesEvals(ida_mem, &ida_njv_tmp);
    if (check_flag(flag, "IDAGetNumJtimesEvals")) { return 1; }
    flag = IDAGetNumLinIters(ida_mem, &ida_nli_tmp);
    if (check_flag(flag, "IDAGetNumLinIters")) { return 1; }
    flag = IDAGetNumErrTestFails(ida_mem, &ida_netf_tmp);
    if (check_flag(flag, "IDAGetNumErrTestFails")) { return 1; }
    flag = IDAGetNumNonlinSolvIters(ida_mem, &ida_nni_tmp);
    if (check_flag(flag, "IDAGetNumNonlinSolvIters")) { return 1; }
    flag = IDAGetNumNonlinSolvConvFails(ida_mem, &ida_ncfn_tmp);
    if (check_flag(flag, "IDAGetNumNonlinSolvConvFails")) { return 1; }

    ida_nst += ida_nst_tmp;
    ida_nre += ida_nre_tmp;
    ida_nsetups += ida_nls_tmp;
    ida_njv += ida_njv_tmp;
    ida_nli += ida_nli_tmp;
    ida_netf += ida_netf_tmp;
    ida_nni += ida_nni_tmp;
    ida_ncfn += ida_ncfn_tmp;
  }

  std::cout << "Final Solver Statistics:\n";
  std::cout << "   Outer PDAEStep steps = " << nst << "\n";
  std::cout << "   Algebraic linear solver setups = " << nlinsetups << "\n";
  std::cout << "   Algebraic nonlinear iterations = " << nniters << "\n";
  std::cout << "   Algebraic nonlinear convergence failures = " << nnfails
            << "\n";
  std::cout << "   Total partition IDA steps = " << ida_nst << "\n";
  std::cout << "   Total partition residual evaluations = " << ida_nre << "\n";
  std::cout << "   Total partition linear solver setups = " << ida_nsetups
            << "\n";
  std::cout << "   Total partition J-times evaluations = " << ida_njv << "\n";
  std::cout << "   Total partition linear iterations = " << ida_nli << "\n";
  std::cout << "   Total partition error test failures = " << ida_netf << "\n";
  std::cout << "   Total partition nonlinear iterations = " << ida_nni << "\n";
  std::cout << "   Total partition nonlinear convergence failures = " << ida_ncfn
            << "\n";

  return 0;
}

static int RunPdae(const UserData& data, const Options& opts, sunrealtype t0,
                   sunrealtype tf, sunrealtype reltol, sunrealtype abstol,
                   int nout, N_Vector uref, SUNContext ctx)
{
  PdaeUserData pdae_data;
  pdae_data.model      = data;
  pdae_data.partitions = opts.partitions;
  pdae_data.coupling_size =
    std::max(static_cast<sunindextype>(1),
             static_cast<sunindextype>(opts.partitions - 1));
  PartitionLayout(data.n, opts.partitions, pdae_data.starts, pdae_data.sizes);

  std::vector<N_Vector> x(opts.partitions, NULL);
  std::vector<N_Vector> z(opts.partitions, NULL);
  std::vector<N_Vector> xp(opts.partitions, NULL);
  std::vector<N_Vector> zp(opts.partitions, NULL);
  std::vector<SUNLinearSolver> partition_ls(opts.partitions, NULL);
  std::vector<void*> ida_mems(opts.partitions, NULL);

  N_Vector w                 = NULL;
  N_Vector wp                = NULL;
  N_Vector y                 = NULL;
  N_Vector yp                = NULL;
  N_Vector u                 = NULL;
  N_Vector alg_template      = NULL;
  N_Vector alg_res           = NULL;
  void* arkode_mem           = NULL;
  SUNNonlinearSolver alg_nls = NULL;
  SUNLinearSolver alg_ls     = NULL;

  std::cout << "\nAdvection pulse problem solved with PDAEStep:\n";
  std::cout << "   PDE                  = u_t = u_x\n";
  std::cout << "   mesh points          = " << data.n << "\n";
  std::cout << "   number of partitions = " << opts.partitions << "\n";
  std::cout << "   dx                   = " << data.dx << "\n";
  std::cout << "   t0                   = " << t0 << "\n";
  std::cout << "   tf                   = " << tf << "\n";
  std::cout << "   pulse start          = " << data.pulse_start << "\n";
  std::cout << "   pulse width          = " << data.pulse_width << "\n";
  std::cout << "   pdae fixed step      = " << opts.pdae_h << "\n";
  std::cout << "   reltol               = " << reltol << "\n";
  std::cout << "   abstol               = " << abstol << "\n\n";

  for (int p = 0; p < opts.partitions; p++)
  {
    x[p] = N_VNew_Serial(pdae_data.sizes[p], ctx);
    if (check_ptr(x[p], "N_VNew_Serial")) { return 1; }
    z[p] = N_VNew_Serial(1, ctx);
    if (check_ptr(z[p], "N_VNew_Serial")) { return 1; }
    xp[p] = N_VClone(x[p]);
    if (check_ptr(xp[p], "N_VClone")) { return 1; }
    zp[p] = N_VClone(z[p]);
    if (check_ptr(zp[p], "N_VClone")) { return 1; }

    N_VConst(ZERO, x[p]);
    N_VConst(ZERO, xp[p]);
    N_VConst((p == opts.partitions - 1) ? Pulse(t0, data) : ZERO, z[p]);
    N_VConst(ZERO, zp[p]);

    sunrealtype* xdata  = N_VGetArrayPointer(x[p]);
    sunrealtype* xpdata = N_VGetArrayPointer(xp[p]);
    sunrealtype* zdata  = N_VGetArrayPointer(z[p]);
    for (sunindextype i = 0; i < pdae_data.sizes[p]; i++)
    {
      sunrealtype uright = (i == pdae_data.sizes[p] - 1) ? zdata[0]
                                                         : xdata[i + 1];
      xpdata[i]          = AdvectionRhsValue(uright, xdata[i], data);
    }
  }

  w = N_VNew_Serial(pdae_data.coupling_size, ctx);
  if (check_ptr(w, "N_VNew_Serial")) { return 1; }
  wp = N_VClone(w);
  if (check_ptr(wp, "N_VClone")) { return 1; }
  N_VConst(ZERO, w);
  N_VConst(ZERO, wp);

  y = PDAEStepManyVector(x.data(), z.data(), w, opts.partitions);
  if (check_ptr(y, "PDAEStepManyVector")) { return 1; }
  yp = PDAEStepManyVector(xp.data(), zp.data(), wp, opts.partitions);
  if (check_ptr(yp, "PDAEStepManyVector")) { return 1; }

  arkode_mem = PDAEStepCreate(PdaeComponentResidualForPartition,
                              PdaeAlgebraicResidual, t0, y, yp, opts.partitions,
                              ctx);
  if (check_ptr(arkode_mem, "PDAEStepCreate")) { return 1; }

  int flag = ARKodeSetUserData(arkode_mem, &pdae_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  flag = PDAEStepSetPartitionJacTimes(arkode_mem,
                                      PdaeComponentJacTimesForPartition);
  if (check_flag(flag, "PDAEStepSetPartitionJacTimes")) { return 1; }

  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  flag = ARKodeSetFixedStep(arkode_mem, opts.pdae_h);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  flag = ARKodeSetMaxNumSteps(arkode_mem, -1);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  flag = PDAEStepSetMaxNonlinIters(arkode_mem, 12);
  if (check_flag(flag, "PDAEStepSetMaxNonlinIters")) { return 1; }

  flag = PDAEStepGetAlgebraicVectorTemplate(arkode_mem, &alg_template);
  if (check_flag(flag, "PDAEStepGetAlgebraicVectorTemplate")) { return 1; }

  alg_nls = SUNNonlinSol_Newton(alg_template, ctx);
  if (check_ptr(alg_nls, "SUNNonlinSol_Newton")) { return 1; }
  flag = PDAEStepSetNonlinearSolver(arkode_mem, alg_nls);
  if (check_flag(flag, "PDAEStepSetNonlinearSolver")) { return 1; }

  alg_ls = SUNLinSol_SPGMR(alg_template, SUN_PREC_NONE, 10, ctx);
  if (check_ptr(alg_ls, "SUNLinSol_SPGMR")) { return 1; }
  flag = PDAEStepSetCouplingJacTimes(arkode_mem, PdaeAlgebraicJacTimes);
  if (check_flag(flag, "PDAEStepSetCouplingJacTimes")) { return 1; }
  flag = PDAEStepSetLinearSolver(arkode_mem, alg_ls, NULL);
  if (check_flag(flag, "PDAEStepSetLinearSolver")) { return 1; }

  for (int p = 0; p < opts.partitions; p++)
  {
    N_Vector partition_template = NULL;
    flag = PDAEStepGetPartitionVectorTemplate(arkode_mem, p, &partition_template);
    if (check_flag(flag, "PDAEStepGetPartitionVectorTemplate")) { return 1; }

    partition_ls[p] = SUNLinSol_SPGMR(partition_template, SUN_PREC_NONE,
                                      pdae_data.sizes[p] + 1, ctx);
    if (check_ptr(partition_ls[p], "SUNLinSol_SPGMR")) { return 1; }

    flag = PDAEStepGetPartitionIntegrator(arkode_mem, p, &ida_mems[p]);
    if (check_flag(flag, "PDAEStepGetPartitionIntegrator")) { return 1; }

    flag = IDASStolerances(ida_mems[p], reltol, abstol);
    if (check_flag(flag, "IDASStolerances")) { return 1; }

    flag = IDASetMaxOrd(ida_mems[p], 1);
    if (check_flag(flag, "IDASetMaxOrd")) { return 1; }

    flag = IDASetLinearSolver(ida_mems[p], partition_ls[p], NULL);
    if (check_flag(flag, "IDASetLinearSolver")) { return 1; }
  }

  u = N_VNew_Serial(data.n, ctx);
  if (check_ptr(u, "N_VNew_Serial")) { return 1; }

  std::ofstream outfile("ark_advection_pulse.txt");
  outfile << "# vars: t x u\n";
  outfile << std::scientific;
  outfile << std::setprecision(std::numeric_limits<sunrealtype>::digits10);

  PrintHeader();
  CopyPdaeSolution(y, pdae_data, u);
  PrintSolution(0, t0, u, data);
  flag = WriteSolution(outfile, t0, u, data);
  if (check_flag(flag, "WriteSolution")) { return 1; }

  sunrealtype tret = t0;
  for (int iout = 1; iout <= nout; iout++)
  {
    sunrealtype tout = tf * static_cast<sunrealtype>(iout) /
                       static_cast<sunrealtype>(nout);

    flag = ARKodeSetStopTime(arkode_mem, tout);
    if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }

    flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
    if (check_flag(flag, "ARKodeEvolve")) { return 1; }

    CopyPdaeSolution(y, pdae_data, u);
    PrintSolution(iout, tret, u, data);
    flag = WriteSolution(outfile, tret, u, data);
    if (check_flag(flag, "WriteSolution")) { return 1; }
  }
  outfile.close();

  alg_res = N_VClone(alg_template);
  if (check_ptr(alg_res, "N_VClone")) { return 1; }
  flag = PdaeAlgebraicResidual(tret, y, w, alg_res, &pdae_data);
  if (check_flag(flag, "PdaeAlgebraicResidual")) { return 1; }

  sunrealtype alg_res_norm = N_VMaxNorm(alg_res);
  std::cout << "\nFinal algebraic residual max norm = " << alg_res_norm << "\n";
  if (uref) { flag = PrintL1Error(u, uref); }

  if (flag == 0 &&
      alg_res_norm > SUNMAX(SUN_RCONST(1.0e-7), SUN_RCONST(100.0) * abstol))
  {
    std::cerr << "ERROR: final algebraic residual norm is too large\n";
    flag = 1;
  }
  else if (flag == 0) { flag = PrintPdaeFinalStats(arkode_mem, ida_mems); }

  ARKodeFree(&arkode_mem);
  if (alg_ls) { SUNLinSolFree(alg_ls); }
  if (alg_nls) { SUNNonlinSolFree(alg_nls); }
  if (alg_res) { N_VDestroy(alg_res); }
  if (u) { N_VDestroy(u); }
  for (int p = 0; p < opts.partitions; p++)
  {
    if (partition_ls[p]) { SUNLinSolFree(partition_ls[p]); }
    N_VDestroy(x[p]);
    N_VDestroy(z[p]);
    N_VDestroy(xp[p]);
    N_VDestroy(zp[p]);
  }
  N_VDestroy(w);
  N_VDestroy(wp);
  N_VDestroy(y);
  N_VDestroy(yp);

  return flag;
}

static int ParseArgs(std::vector<std::string>& args, UserData& data,
                     sunrealtype& tf, sunrealtype& reltol, sunrealtype& abstol,
                     int& nout, Options& opts)
{
  find_arg(args, "--n", data.n);
  find_arg(args, "--tf", tf);
  find_arg(args, "--rtol", reltol);
  find_arg(args, "--atol", abstol);
  find_arg(args, "--nout", nout);
  find_arg(args, "--pulse-start", data.pulse_start);
  find_arg(args, "--pulse-width", data.pulse_width);
  find_arg(args, "--solver", opts.solver);
  find_arg(args, "--partitions", opts.partitions);
  find_arg(args, "--ida-max-order", opts.ida_max_order);
  find_arg(args, "--pdae-h", opts.pdae_h);

  if (!args.empty())
  {
    std::cerr << "ERROR: unknown argument '" << args.front() << "'\n";
    return 1;
  }

  if (data.n < 2)
  {
    std::cerr << "ERROR: --n must be at least 2\n";
    return 1;
  }

  if (tf <= ZERO)
  {
    std::cerr << "ERROR: --tf must be positive\n";
    return 1;
  }

  if (nout < 1)
  {
    std::cerr << "ERROR: --nout must be positive\n";
    return 1;
  }

  if (data.pulse_width <= ZERO)
  {
    std::cerr << "ERROR: --pulse-width must be positive\n";
    return 1;
  }

  if (opts.solver != "ida" && opts.solver != "pdae")
  {
    std::cerr << "ERROR: --solver must be 'ida' or 'pdae'\n";
    return 1;
  }

  if (opts.partitions < 1)
  {
    std::cerr << "ERROR: --partitions must be positive\n";
    return 1;
  }

  if (opts.partitions > data.n)
  {
    std::cerr << "ERROR: --partitions must not exceed --n\n";
    return 1;
  }

  if (opts.ida_max_order < 0)
  {
    std::cerr << "ERROR: --ida-max-order must be nonnegative\n";
    return 1;
  }

  if (opts.pdae_h <= ZERO)
  {
    std::cerr << "ERROR: --pdae-h must be positive\n";
    return 1;
  }

  data.dx = ONE / static_cast<sunrealtype>(data.n);

  return 0;
}
