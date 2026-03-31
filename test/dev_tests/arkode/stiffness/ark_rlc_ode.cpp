#include <cstdlib>
#include <cstring>
#include <fstream>

// Problem header
#include <problems/rlc_ode.hpp>

// Test utilities
#include <utilities/parse_options.hpp>
#include <utilities/compute_eigenvalues.hpp>
#include <utilities/compute_eigenvalues_2.hpp>

// SUNDIALS headers
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "arkode/arkode.h"

using namespace std;
using namespace problems::rlc_ode;

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Default problem parameters
  sunrealtype t0             = 0.0;
  sunrealtype tf             = 15000.0;
  sunrealtype dt_out         = 0.1;
  string output_prefix       = "rlc_ode";
  bool stiffness_detection   = false;
  sunrealtype voltage_source = 12.0;

  vector<OptionSpec> options;

  options.push_back({"-v", "--voltage-source", true,
                     "Set beta parameter (default: 12.0)", [&](const char* v)
                     { return setReal(voltage_source, v); }});

  options.push_back({"-s", "--stiffness-detection", false,
                     "Enable stiffness detection (default: false)",
                     [&](const char*)
                     { return setBoolTrue(stiffness_detection); }});

  options.push_back({"-p", "--prefix", true,
                     "Output file name prefix (default: rlc_ode)",
                     [&](const char* v) { return setString(output_prefix, v); }});

  if (parseArguments(argc, argv, options)) { return 1; }

  // Error flag
  int ierr = 0;

  // Create problem instance
  ODEProblem problem;
  problem.setVoltageSource(voltage_source);
  int NEQ = problem.getNumEquations();

  // Create SUNDIALS vector for initial conditions
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  if (y == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    return 1;
  }

  // Set initial condition
  N_VConst(0.0, y);

  // Create ARKODE memory structure
  void* arkode_mem = ARKStepCreate(ODEProblem::rhsWrapper, nullptr, t0, y,
                                   sunctx);
  if (arkode_mem == nullptr)
  {
    cerr << "Error creating ARKODE memory" << endl;
    return 1;
  }

  if (stiffness_detection)
  {
    ierr = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE",
                               "ARKODE_FEHLBERG_SHAMPINE_HIEBERT_6_4_5");
    if (ierr)
    {
      cerr << "Error setting method table" << endl;
      return 1;
    }
  }

  // Set user data (pointer to problem instance)
  ARKodeSetUserData(arkode_mem, &problem);

  // Set tolerances
  sunrealtype reltol = 1e-8;
  sunrealtype abstol = 1e-10;
  ARKodeSStolerances(arkode_mem, reltol, abstol);

  SUNMatrix M        = SUNDenseMatrix(NEQ, NEQ, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, M, sunctx);
  ARKodeSetMassLinearSolver(arkode_mem, LS, M, false);
  ARKodeSetMassFn(arkode_mem, ODEProblem::massDenseWrapper);

  SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  sunrealtype* Adata = SUNDenseMatrix_Data(A);

  sunrealtype R1 = problem.getR1();
  sunrealtype R2 = problem.getR2();
  sunrealtype R3 = problem.getR3();
  sunrealtype R4 = problem.getR4();
  sunrealtype L2 = problem.getL2();
  sunrealtype L4 = problem.getL4();
  sunrealtype C3 = problem.getC3();
  sunrealtype C4 = problem.getC4();
  sunrealtype Vs = problem.getVs();

  // Column 0 - d/di2
  Adata[0] = (-(R1 + R2) - R3) / L2;
  Adata[1] = 1.0 / C3;
  Adata[2] = R3 / L4;
  Adata[3] = 0.0;
  // Column 1 - d/dvC3
  Adata[4] = -1.0 / L2;
  Adata[5] = 0.0;
  Adata[6] = 1.0 / L4;
  Adata[7] = 0.0;
  // Column 2 - d/di4
  Adata[8]  = R3 / L2;
  Adata[9]  = -1.0 / C3;
  Adata[10] = (-R3 - R4) / L4;
  Adata[11] = 1.0 / C4;
  // Column 3 - d/dvC4
  Adata[12] = 0.0;
  Adata[13] = 0.0;
  Adata[14] = -1.0 / L4;
  Adata[15] = 0.0;

  cout << "A matrix:\n";
  SUNDenseMatrix_Print(A, stdout);
  cout << "\n";

  auto eigs = computeEigenvalues(Adata); // will modify Adata
  ofstream eigfile(output_prefix + "_eig.txt");
  eigs.print();
  eigs.print(eigfile, true);
  cout << eigs.stiffness_ratio << "\n";

  cout << "A matrix:\n";
  SUNDenseMatrix_Print(A, stdout);
  cout << "\n";

  problem.computeMassDense(M);
  problem.computeJacDense(A);
  SUNMatrix J = form_jacobian_from_mass(M, A);
  cout << "A matrix:\n";
  SUNDenseMatrix_Print(J, stdout);


  // Time integration loop
  sunrealtype t    = t0;
  sunrealtype tout = t0 + dt_out;
  int steps        = 0;

  // Open file for writing data
  ofstream datafile(output_prefix + "_data.txt");
  problem.printHeader({"i2", "vC3", "i4", "vC4"}, datafile);
  problem.printState(t, y, datafile);

  problem.printHeader({"i2", "vC3", "i4", "vC4"});
  problem.printState(t, y);

  while (t < tf)
  {
    int flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_ONE_STEP);
    if (flag < 0)
    {
      cerr << "ARKODE error, flag = " << flag << endl;
      break;
    }
    steps++;

    problem.printState(t, y, datafile);
    if (steps % 1000 == 0) { problem.printState(t, y); }

    tout += dt_out;
    if (tout > tf) tout = tf;
  }

  datafile.close();

  // Print solver statistics
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // Clean up
  ARKodeFree(&arkode_mem);
  N_VDestroy(y);

  return 0;
}
