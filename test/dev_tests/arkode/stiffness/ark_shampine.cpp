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
void printUsage(const char* progname)
{
  cout << "Usage: " << progname << " [options]\n\n";
  cout << "Options:\n";
  cout << "  -b, --beta <value>     Set beta parameter (default: 5.0)\n";
  cout << "  -g, --gamma <value>    Set gamma parameter (default: 2.0)\n";
  cout << "  -o, --output <file>    Output file name (default: output to "
          "console)\n";
  cout << "  -h, --help             Print this help message\n";
  cout << "\nExample:\n";
  cout << "  " << progname << " --beta 3.5 --gamma 1.5 --output results.txt\n";
  cout << "  " << progname << " -b 10 -g 5 -o data.dat\n";
}

// Function to parse command line arguments
bool parseArguments(int argc, char* argv[], sunrealtype& beta,
                    sunrealtype& gamma, string& outputFile)
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

  string output_file = "data.txt";

  // Parse command line arguments
  if (!parseArguments(argc, argv, beta, gamma, output_file))
  {
    return (argc > 1 &&
            (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
             ? 0
             : 1;
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

  // Create ARKODE memory structure
  void* arkode_mem = ARKStepCreate(nullptr, ODEProblem::rhsWrapper, t0, y,
                                   sunctx);
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

  // Create dense matrix and linear solver
  SUNMatrix A = SUNDenseMatrix(problem.getNumEquations(),
                               problem.getNumEquations(), sunctx);
  if (A == nullptr)
  {
    cerr << "Error creating matrix" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
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

  // Open file for writing data
  ofstream datafile(output_file);
  datafile << setprecision(17) << scientific;
  datafile << "# t, y1, y2, y3\n";

  // Output initial condition
  sunrealtype* ydata = N_VGetArrayPointer(y);
  datafile << setw(22) << t0 << setw(25) << ydata[0] << setw(25) << ydata[1]
           << setw(25) << ydata[2] << endl;

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
             << setw(25) << ydata[2] << endl;

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
