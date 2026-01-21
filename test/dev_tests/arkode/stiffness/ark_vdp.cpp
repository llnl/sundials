#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// Problem header
#include <problems/van_der_pol.hpp>

// SUNDIALS headers
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "arkode/arkode.h"

using namespace std;
using namespace problems::van_der_pol;

// Function to print usage information
void printUsage(const char* progname)
{
  cout << "Usage: " << progname << " [options]\n\n";
  cout << "Solves the Van der Pol equation:\n";
  cout << "  y'' - mu*(1 - y^2)*y' + y = 0\n\n";
  cout << "Options:\n";
  cout << "  -m, --mu <value>       Set mu parameter (default: 100.0)\n";
  cout << "  -e, --explicit         Use ERK method (default: DIRK)\n";
  cout << "  -o, --output <file>    Output file name (default: data.txt)\n";
  cout << "  -h, --help             Print this help message\n";
}

// Function to parse command line arguments
bool parseArguments(int argc, char* argv[], sunrealtype& mu, string& outputFile,
                    bool& use_explicit)
{
  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mu") == 0)
    {
      if (i + 1 < argc) { mu = atof(argv[++i]); }
      else
      {
        cerr << "Error: --mu requires a value\n";
        return false;
      }
    }
    else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--explicit") == 0)
    {
      use_explicit = true;
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
  sunrealtype mu     = 100.0;
  sunrealtype t0     = 0.0;
  sunrealtype tf     = 1000.0;
  sunrealtype dt_out = 0.1;

  string output_file = "data.txt";
  bool use_explicit  = false;

  // Parse command line arguments
  if (!parseArguments(argc, argv, mu, output_file, use_explicit))
  {
    return (argc > 1 &&
            (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
             ? 0
             : 1;
  }

  // Create problem instance
  ODEProblem problem(mu);

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
  datafile << "# t, y1, y2\n";

  // Output initial condition
  sunrealtype* ydata = N_VGetArrayPointer(y);
  datafile << setw(22) << t0 << setw(25) << ydata[0] << setw(25) << ydata[1]
           << endl;

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

    // Print solution
    datafile << setw(22) << t << setw(25) << ydata[0] << setw(25) << ydata[1]
             << endl;

    tout += dt_out;
    if (tout > tf) tout = tf;
  }

  datafile.close();

  // Print solver statistics
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // Clean up
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);

  return 0;
}
