#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// Problem header
#include <problems/shampine.hpp>

// SUNDIALS headers
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "arkode/arkode.h"

using namespace std;
using namespace problems::shampine;

// Function to print usage information
static void printUsage(const char* progname)
{
  cout << "Usage: " << progname << " [options]\n\n";
  cout << "Options:\n";
  cout << "  -b, --beta <value>     Set beta parameter (default: 5.0)\n";
  cout << "  -g, --gamma <value>    Set gamma parameter (default: 2.0)\n";
  cout << "  -e, --explicit         Use ERK method (default: DIRK)\n";
  cout << "  -c, --comparison       Use comparison method (default: false)\n";
  cout << "  -o, --output <file>    Output file name (default: data.txt)\n";
  cout << "  -h, --help             Print this help message\n";
}

// Function to parse command line arguments
static bool parseArguments(int argc, char* argv[], sunrealtype& beta,
                           sunrealtype& gamma, string& outputFile,
                           bool& use_explicit, bool& use_comparison)
{
  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--beta") == 0)
    {
      if (i + 1 < argc) { beta = atof(argv[++i]); }
      else
      {
        cerr << "Error: --beta requires a value\n";
        return false;
      }
    }
    else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--gamma") == 0)
    {
      if (i + 1 < argc) { gamma = atof(argv[++i]); }
      else
      {
        cerr << "Error: --gamma requires a value\n";
        return false;
      }
    }
    else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--explicit") == 0)
    {
      use_explicit = true;
    }
    else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--comparison") == 0)
    {
      use_comparison = true;
    }
    else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0)
    {
      if (i + 1 < argc) { outputFile = argv[++i]; }
      else
      {
        cerr << "Error: --output requires a filename\n";
        return false;
      }
    }
    else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
    {
      printUsage(argv[0]);
      return false;
    }
    else
    {
      cerr << "Error: Unknown option '" << argv[i] << "'\n\n";
      printUsage(argv[0]);
      return false;
    }
  }
  return true;
}

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Default problem parameters
  sunrealtype beta   = 10.0;
  sunrealtype gamma  = 10.0;
  sunrealtype t0     = 0.0;
  sunrealtype tf     = 1.0;
  sunrealtype dt_out = 0.1;

  string output_file  = "data.txt";
  bool use_explicit   = false;
  bool use_comparison = false;

  int ierr = 0;

  // Parse command line arguments
  if (!parseArguments(argc, argv, beta, gamma, output_file, use_explicit,
                      use_comparison))
  {
    return (argc > 1 &&
            (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
             ? 0
             : 1;
  }

  if (use_comparison && !use_explicit)
  {
    cerr << "Error use comparison requires use explicit" << endl;
    return 1;
  }

  // Create problem instance
  ODEProblem problem(beta, gamma);

  // Create SUNDIALS vector for initial conditions
  N_Vector y = N_VNew_Serial(problem.getNumEquations(), sunctx);
  if (y == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    return 1;
  }

  // Set initial conditions
  problem.setInitialConditions(y);

  // Create ARKODE memory structure based on method choice
  void* arkode_mem = nullptr;
  if (use_explicit)
  {
    cout << "Using ERK method" << endl;
    arkode_mem = ARKStepCreate(ODEProblem::rhsWrapper, nullptr, t0, y, sunctx);
  }
  else
  {
    cout << "Using DIRK method" << endl;
    arkode_mem = ARKStepCreate(nullptr, ODEProblem::rhsWrapper, t0, y, sunctx);
  }

  if (arkode_mem == nullptr)
  {
    cerr << "Error creating ARKODE memory" << endl;
    N_VDestroy(y);
    return 1;
  }

  if (use_comparison)
  {
    ierr = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", "ARKODE_FEHLBERG_SHAMPINE_HIEBERT_6_4_5");
    if (ierr)
    {
      cerr << "Error setting method table" << endl;
      return 1;
    }
  }

  // Set user data (pointer to problem instance)
  ARKodeSetUserData(arkode_mem, &problem);

  // Set tolerances
  sunrealtype reltol = 1e-6;
  sunrealtype abstol = 1e-8;
  ARKodeSStolerances(arkode_mem, reltol, abstol);

  // Only set up linear solver and Jacobian for implicit method
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (!use_explicit)
  {
    // Create dense matrix and linear solver
    A = SUNDenseMatrix(problem.getNumEquations(), problem.getNumEquations(),
                       sunctx);
    if (A == nullptr)
    {
      cerr << "Error creating matrix" << endl;
      ARKodeFree(&arkode_mem);
      N_VDestroy(y);
      return 1;
    }

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (LS == nullptr)
    {
      cerr << "Error creating linear solver" << endl;
      SUNMatDestroy(A);
      ARKodeFree(&arkode_mem);
      N_VDestroy(y);
      return 1;
    }

    // Attach linear solver
    ARKodeSetLinearSolver(arkode_mem, LS, A);

    // Set Jacobian function
    ARKodeSetJacFn(arkode_mem, ODEProblem::jacWrapper);
  }

  // Open file for writing data
  ofstream datafile(output_file);
  datafile << setprecision(17) << scientific;
  datafile << "# t, y1, y2, y3, est err1, est err2, est err3, lambda1, "
              "lambda2, lambda3, stiffness ratio\n";

  // Initial output
  sunrealtype* ydata = N_VGetArrayPointer(y);

  complex<sunrealtype> lambda1, lambda2, lambda3;
  problem.computeEigenvalues(lambda1, lambda2, lambda3);
  auto stiffness_ratio = problem.computeStiffnessRatio();

  datafile
  << setw(26) << t0
  << setw(26) << ydata[0]
  << setw(26) << ydata[1]
  << setw(26) << ydata[2]
  << setw(26) << 0.0
  << setw(26) << 0.0
  << setw(26) << 0.0
  << setw(26) << 0.0
  << setw(26) << 0.0
  << setw(26) << 0.0
  << setw(26) << lambda1.real() << showpos << lambda1.imag() << noshowpos << "j"
  << setw(26) << lambda2.real() << showpos << lambda2.imag() << noshowpos << "j"
  << setw(26) << lambda3.real() << showpos << lambda3.imag() << noshowpos << "j"
  << setw(26) << stiffness_ratio
  << endl;

  // Vector for local error estimate
  N_Vector loc_err_est = N_VClone(y);
  if (loc_err_est == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  N_Vector companion_loc_err_est = N_VClone(y);
  if (companion_loc_err_est == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    N_VDestroy(loc_err_est);
    return 1;
  }

  N_Vector err_weights = N_VClone(y);
  if (err_weights == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    N_VDestroy(loc_err_est);
    N_VDestroy(companion_loc_err_est);
    return 1;
  }

  // Time integration loop
  sunrealtype t    = t0;
  sunrealtype tout = t0 + dt_out;

  while (t < tf)
  {
    int flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_ONE_STEP);

    if (flag < 0)
    {
      cerr << "ARKODE error, flag = " << flag << endl;
      break;
    }

    // Step output
    flag = ARKodeGetEstLocalErrors(arkode_mem, loc_err_est);
    if (flag < 0)
    {
      cerr << "ARKODE error, flag = " << flag << endl;
      break;
    }
    sunrealtype* lee_data = N_VGetArrayPointer(loc_err_est);

    flag = ARKodeGetEstLocalErrors2(arkode_mem, companion_loc_err_est);
    if (flag < 0)
    {
      cerr << "ARKODE error, flag = " << flag << endl;
      break;
    }
    sunrealtype* clee_data = N_VGetArrayPointer(companion_loc_err_est);

    flag = ARKodeGetErrWeights(arkode_mem, err_weights);
    if (flag < 0)
    {
      cerr << "ARKODE error, flag = " << flag << endl;
      break;
    }
    sunrealtype* weights_data = N_VGetArrayPointer(err_weights);

    problem.computeEigenvalues(lambda1, lambda2, lambda3);
    stiffness_ratio = problem.computeStiffnessRatio();

    datafile
    << setw(26) << t
    << setw(26) << ydata[0]
    << setw(26) << ydata[1]
    << setw(26) << ydata[2]
    << setw(26) << lee_data[0] * weights_data[0]
    << setw(26) << lee_data[1] * weights_data[1]
    << setw(26) << lee_data[2] * weights_data[2]
    << setw(26) << clee_data[0] * weights_data[0]
    << setw(26) << clee_data[1] * weights_data[1]
    << setw(26) << clee_data[2] * weights_data[2]
    << setw(26) << lambda1.real() << showpos << lambda1.imag() << noshowpos << "j"
    << setw(26) << lambda2.real() << showpos << lambda2.imag() << noshowpos << "j"
    << setw(26) << lambda3.real() << showpos << lambda3.imag() << noshowpos << "j"
    << setw(26) << stiffness_ratio
    << endl;

    tout += dt_out;
    if (tout > tf) tout = tf;
  }

  datafile.close();

  // Print solver statistics
  cout << "Eigenvalue 1                  = " << lambda1.real() << showpos
  << lambda1.imag() << noshowpos << "i" << " (mag = " << abs(lambda1) << ")" << endl;
  cout << "Eigenvalue 2                  = " << lambda2.real() << showpos
       << lambda2.imag() << noshowpos << "i"  << " (mag = " << abs(lambda2) << ")" << endl;
  cout << "Eigenvalue 3                  = " << lambda3.real() << endl;
  cout << "Stiffness ratio               = " << stiffness_ratio << endl;
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // Clean up
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(loc_err_est);
  N_VDestroy(companion_loc_err_est);
  N_VDestroy(err_weights);

  return 0;
}
