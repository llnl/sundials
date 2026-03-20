/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * Header file for ARKODE LSRKStep Sod shock tube example, see
 * ark_sod_lsrk.cpp for more details.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_lsrkstep.h"
#include "nvector/nvector_manyvector.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_core.hpp"

// Macros for problem constants
#define rhoL   SUN_RCONST(1.0)
#define rhoR   SUN_RCONST(0.125)
#define pL     SUN_RCONST(1.0)
#define pR     SUN_RCONST(0.1)
#define uL     SUN_RCONST(0.0)
#define uR     SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define FOURTH SUN_RCONST(0.25)

#define NSPECIES 5
#define STSIZE   6

#define WIDTH (10 + std::numeric_limits<sunrealtype>::digits10)

// -----------------------------------------------------------------------------
// Problem options
// -----------------------------------------------------------------------------

class ARKODEParameters
{
public:
  // Integration method (ARKODE_LSRK_SSP_S_2, ARKODE_LSRK_SSP_S_3, ARKODE_LSRK_SSP_10_4,
  // or any valid ARKODE_ERKTableID for ERK methods)
  std::string integrator;

  // Method stages (0 => to use the default; ignored if using ARKODE_LSRK_SSP_10_4 or
  // an ERK method)
  int stages;

  // Relative and absolute tolerances
  sunrealtype rtol;
  sunrealtype atol;

  // Step size selection (ZERO = adaptive steps)
  sunrealtype fixed_h;

  // Maximum number of time steps between outputs
  int maxsteps;

  // Output-related information
  int output;         // 0 = none, 1 = stats, 2 = disk, 3 = disk with tstop
  int nout;           // number of output times
  std::ofstream uout; // output file stream

  // constructor (with default values)
  ARKODEParameters()
    : integrator("ARKODE_LSRK_SSP_10_4"),
      stages(0),
      rtol(SUN_RCONST(1.e-4)),
      atol(SUN_RCONST(1.e-11)),
      fixed_h(ZERO),
      maxsteps(10000),
      output(1),
      nout(10){};

}; // end ARKODEParameters

// -----------------------------------------------------------------------------
// Problem parameters
// -----------------------------------------------------------------------------

// user data class
class EulerData
{
public:
  ///// domain related data /////
  long int nx;    // global number of x grid points
  sunrealtype t0; // time domain extents
  sunrealtype tf;
  sunrealtype xl; // spatial domain extents
  sunrealtype xr;
  sunrealtype dx; // spatial mesh spacing

  ///// problem-defining data /////
  sunrealtype gamma; // ratio of specific heat capacities, cp/cv

  ///// reusable arrays for WENO flux calculations /////
  sunrealtype* flux;
  sunrealtype w1d[STSIZE][NSPECIES];

  ///// class operations /////

  // constructor
  EulerData()
    : nx(512),
      t0(ZERO),
      tf(SUN_RCONST(0.1)),
      xl(ZERO),
      xr(ONE),
      dx(ZERO),
      gamma(SUN_RCONST(1.4)),
      flux(nullptr){};

  // manual destructor
  void FreeData()
  {
    delete[] flux;
    flux = nullptr;
  };

  // destructor
  ~EulerData() { this->FreeData(); };

  // Utility routine to pack 1-dimensional data for *interior only* data;
  // e.g., in the x-direction given a location (i), we return values at
  // the 6 nodal values closest to the (i-1/2) face along the x-direction,
  // {w(i-3), w(i-2), w(i-1), w(i), w(i+1), w(i+2)}.
  inline void pack1D(const sunrealtype* rho, const sunrealtype* mx,
                     const sunrealtype* my, const sunrealtype* mz,
                     const sunrealtype* et, const long int& i)
  {
    for (int l = 0; l < STSIZE; l++) { this->w1d[l][0] = rho[i - 3 + l]; }
    for (int l = 0; l < STSIZE; l++) { this->w1d[l][1] = mx[i - 3 + l]; }
    for (int l = 0; l < STSIZE; l++) { this->w1d[l][2] = my[i - 3 + l]; }
    for (int l = 0; l < STSIZE; l++) { this->w1d[l][3] = mz[i - 3 + l]; }
    for (int l = 0; l < STSIZE; l++) { this->w1d[l][4] = et[i - 3 + l]; }
  }

  // Utility routine to pack 1-dimensional data for locations near the
  // boundary; like the routine above this packs the 6 closest
  // entries aligned with, e.g., the (i-1/2) face, but now some entries
  // are set to satisfy homogeneous Neumann boundary conditions.
  inline void pack1D_bdry(const sunrealtype* rho, const sunrealtype* mx,
                          const sunrealtype* my, const sunrealtype* mz,
                          const sunrealtype* et, const long int& i)
  {
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l][0] = (i < (3 - l)) ? rho[2 - (i + l)] : rho[i - 3 + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l][1] = (i < (3 - l)) ? mx[2 - (i + l)] : mx[i - 3 + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l][2] = (i < (3 - l)) ? my[2 - (i + l)] : my[i - 3 + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l][3] = (i < (3 - l)) ? mz[2 - (i + l)] : mz[i - 3 + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l][4] = (i < (3 - l)) ? et[2 - (i + l)] : et[i - 3 + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l + 3][0] = (i > (nx - l - 1)) ? rho[i + l - 3] : rho[i + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l + 3][1] = (i > (nx - l - 1)) ? mx[i + l - 3] : mx[i + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l + 3][2] = (i > (nx - l - 1)) ? my[i + l - 3] : my[i + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l + 3][3] = (i > (nx - l - 1)) ? mz[i + l - 3] : mz[i + l];
    }
    for (int l = 0; l < 3; l++)
    {
      this->w1d[l + 3][4] = (i > (nx - l - 1)) ? et[i + l - 3] : et[i + l];
    }
  }

  // Equation of state -- compute and return pressure,
  //    p = (gamma-1)*(e - rho/2*(vx^2+vy^2+vz^2)), or equivalently
  //    p = (gamma-1)*(e - (mx^2+my^2+mz^2)/(2*rho))
  inline sunrealtype eos(const sunrealtype& rho, const sunrealtype& mx,
                         const sunrealtype& my, const sunrealtype& mz,
                         const sunrealtype& et) const
  {
    return ((gamma - ONE) * (et - (mx * mx + my * my + mz * mz) * HALF / rho));
  }

  // Equation of state inverse -- compute and return energy,
  //    e_t = p/(gamma-1) + rho/2*(vx^2+vy^2+vz^2), or equivalently
  //    e_t = p/(gamma-1) + (mx^2+my^2+mz^2)/(2*rho)
  inline sunrealtype eos_inv(const sunrealtype& rho, const sunrealtype& mx,
                             const sunrealtype& my, const sunrealtype& mz,
                             const sunrealtype& pr) const
  {
    return (pr / (gamma - ONE) + (mx * mx + my * my + mz * mz) * HALF / rho);
  }

}; // end EulerData;

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// ODE right hand side (RHS) functions
int frhs(sunrealtype t, N_Vector y, N_Vector f, void* user_data);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// WENO flux calculation helper function
void face_flux(sunrealtype (&w1d)[6][NSPECIES], sunrealtype* f_face,
               const EulerData& udata);

// Compute the initial condition
int SetIC(N_Vector y, EulerData& udata);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Check function return flag
static int check_flag(int flag, const std::string funcname)
{
  if (flag < 0)
  {
    std::cerr << "ERROR: " << funcname << " returned " << flag << std::endl;
    return 1;
  }
  return 0;
}

// Check if a function returned a NULL pointer
static int check_ptr(void* ptr, const std::string funcname)
{
  if (ptr) { return 0; }
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

// Print command line options
static void InputHelp()
{
  std::cout << std::endl;
  std::cout << "Command line options:" << std::endl;
  std::cout << "  --integrator <str> : method (ARKODE_LSRK_SSP_S_2, "
               "ARKODE_LSRK_SSP_S_3, "
               "ARKODE_LSRK_SSP_10_4, or any valid ARKODE_ERKTableID)\n";
  std::cout << "  --stages <int>     : number of stages (ignored for "
               "ARKODE_LSRK_SSP_10_4 and ERK)\n";
  std::cout << "  --tf <real>        : final time\n";
  std::cout << "  --xl <real>        : domain lower boundary\n";
  std::cout << "  --xr <real>        : domain upper boundary\n";
  std::cout << "  --gamma <real>     : ideal gas constant\n";
  std::cout << "  --nx <int>         : number of mesh points\n";
  std::cout << "  --rtol <real>      : relative tolerance\n";
  std::cout << "  --atol <real>      : absolute tolerance\n";
  std::cout << "  --fixed_h <real>   : fixed step size\n";
  std::cout << "  --maxsteps <int>   : max steps between outputs\n";
  std::cout << "  --output <int>     : output level\n";
  std::cout << "  --nout <int>       : number of outputs\n";
  std::cout << "  --help             : print options and exit\n";
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     sunrealtype& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    dest = stof(*(it + 1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    dest = stod(*(it + 1));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    dest = stold(*(it + 1));
#endif
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     long int& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     int& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     std::string& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = *(it + 1);
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     bool& dest, bool store = true)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}

static int ReadInputs(std::vector<std::string>& args, EulerData& udata,
                      ARKODEParameters& uopts, SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem parameters
  find_arg(args, "--gamma", udata.gamma);
  find_arg(args, "--tf", udata.tf);
  find_arg(args, "--xl", udata.xl);
  find_arg(args, "--xr", udata.xr);
  find_arg(args, "--nx", udata.nx);

  // Integrator options
  find_arg(args, "--integrator", uopts.integrator);
  find_arg(args, "--stages", uopts.stages);
  find_arg(args, "--rtol", uopts.rtol);
  find_arg(args, "--atol", uopts.atol);
  find_arg(args, "--fixed_h", uopts.fixed_h);
  find_arg(args, "--maxsteps", uopts.maxsteps);
  find_arg(args, "--output", uopts.output);
  find_arg(args, "--nout", uopts.nout);

  // Recompute mesh spacing and [re]allocate flux array
  udata.dx = (udata.xr - udata.xl) / ((sunrealtype)udata.nx);
  if (udata.flux) { delete[] udata.flux; }
  udata.flux = new sunrealtype[NSPECIES * (udata.nx + 1)];

  if (uopts.stages < 0)
  {
    std::cerr << "ERROR: Invalid number of stages" << std::endl;
    return -1;
  }

  return 0;
}

// Print user data
static int PrintSetup(EulerData& udata, ARKODEParameters& uopts)
{
  std::cout << std::endl;
  std::cout << "Problem parameters and options:" << std::endl;
  std::cout << " --------------------------------- " << std::endl;
  std::cout << "  gamma      = " << udata.gamma << std::endl;
  std::cout << " --------------------------------- " << std::endl;
  std::cout << "  tf         = " << udata.tf << std::endl;
  std::cout << "  xl         = " << udata.xl << std::endl;
  std::cout << "  xr         = " << udata.xr << std::endl;
  std::cout << "  nx         = " << udata.nx << std::endl;
  std::cout << "  dx         = " << udata.dx << std::endl;
  std::cout << " --------------------------------- " << std::endl;
  std::cout << "  integrator = " << uopts.integrator << std::endl;
  if (uopts.stages > 0)
  {
    std::cout << "  stages     = " << uopts.stages << std::endl;
  }
  std::cout << "  rtol       = " << uopts.rtol << std::endl;
  std::cout << "  atol       = " << uopts.atol << std::endl;
  std::cout << "  fixed h    = " << uopts.fixed_h << std::endl;
  std::cout << " --------------------------------- " << std::endl;
  std::cout << "  output     = " << uopts.output << std::endl;
  std::cout << " --------------------------------- " << std::endl;
  std::cout << std::endl;

  return 0;
}

// Initialize output
static int OpenOutput(EulerData& udata, ARKODEParameters& uopts)
{
  // Header for status output
  if (uopts.output)
  {
    std::cout << std::scientific;
    std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
    std::cout << "    t        "
              << " ||rho||     "
              << " ||mx||      "
              << " ||my||      "
              << " ||mz||      "
              << " ||et||" << std::endl;
    std::cout
      << " -----------------------------------------------------------------"
         "---------"
      << std::endl;
  }

  // Open output stream and output problem information
  if (uopts.output >= 2)
  {
    // Open output stream
    std::stringstream fname;
    fname << "sod.out";
    uopts.uout.open(fname.str());

    uopts.uout << std::scientific;
    uopts.uout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
    uopts.uout << "# title Sod Shock Tube" << std::endl;
    uopts.uout << "# nvar 5" << std::endl;
    uopts.uout << "# vars rho mx my mz et" << std::endl;
    uopts.uout << "# nt " << uopts.nout + 1 << std::endl;
    uopts.uout << "# nx " << udata.nx << std::endl;
    uopts.uout << "# xl " << udata.xl << std::endl;
    uopts.uout << "# xr " << udata.xr << std::endl;
  }

  return 0;
}

// Write output
static int WriteOutput(sunrealtype t, N_Vector y, EulerData& udata,
                       ARKODEParameters& uopts)
{
  if (uopts.output)
  {
    // Compute rms norm of the state
    N_Vector rho       = N_VGetSubvector_ManyVector(y, 0);
    N_Vector mx        = N_VGetSubvector_ManyVector(y, 1);
    N_Vector my        = N_VGetSubvector_ManyVector(y, 2);
    N_Vector mz        = N_VGetSubvector_ManyVector(y, 3);
    N_Vector et        = N_VGetSubvector_ManyVector(y, 4);
    sunrealtype rhorms = sqrt(N_VDotProd(rho, rho) / (sunrealtype)udata.nx);
    sunrealtype mxrms  = sqrt(N_VDotProd(mx, mx) / (sunrealtype)udata.nx);
    sunrealtype myrms  = sqrt(N_VDotProd(my, my) / (sunrealtype)udata.nx);
    sunrealtype mzrms  = sqrt(N_VDotProd(mz, mz) / (sunrealtype)udata.nx);
    sunrealtype etrms  = sqrt(N_VDotProd(et, et) / (sunrealtype)udata.nx);
    std::cout << std::setprecision(2) << "  " << t << std::setprecision(5)
              << "  " << rhorms << "  " << mxrms << "  " << myrms << "  "
              << mzrms << "  " << etrms << std::endl;

    // Write solution to disk
    if (uopts.output >= 2)
    {
      sunrealtype* rhodata = N_VGetArrayPointer(rho);
      if (check_ptr(rhodata, "N_VGetArrayPointer")) { return -1; }
      sunrealtype* mxdata = N_VGetArrayPointer(mx);
      if (check_ptr(mxdata, "N_VGetArrayPointer")) { return -1; }
      sunrealtype* mydata = N_VGetArrayPointer(my);
      if (check_ptr(mydata, "N_VGetArrayPointer")) { return -1; }
      sunrealtype* mzdata = N_VGetArrayPointer(mz);
      if (check_ptr(mzdata, "N_VGetArrayPointer")) { return -1; }
      sunrealtype* etdata = N_VGetArrayPointer(et);
      if (check_ptr(etdata, "N_VGetArrayPointer")) { return -1; }

      uopts.uout << t;
      for (sunindextype i = 0; i < udata.nx; i++)
      {
        uopts.uout << std::setw(WIDTH) << rhodata[i];
        uopts.uout << std::setw(WIDTH) << mxdata[i];
        uopts.uout << std::setw(WIDTH) << mydata[i];
        uopts.uout << std::setw(WIDTH) << mzdata[i];
        uopts.uout << std::setw(WIDTH) << etdata[i];
      }
      uopts.uout << std::endl;
    }
  }

  return 0;
}

// Finalize output
static int CloseOutput(ARKODEParameters& uopts)
{
  // Footer for status output
  if (uopts.output)
  {
    std::cout
      << " -----------------------------------------------------------------"
         "---------"
      << std::endl;
    std::cout << std::endl;
  }

  // Close output streams
  if (uopts.output >= 2) { uopts.uout.close(); }

  return 0;
}

//---- end of file ----
