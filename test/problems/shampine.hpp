#ifndef SHAMPINE_HPP_
#define SHAMPINE_HPP_

#include <complex>
#include <limits>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace problems {
namespace shampine {

class ODEProblem
{
private:
  sunrealtype beta;
  sunrealtype gamma;
  static const int NEQ = 3; // Number of equations

public:
  // Constructor
  ODEProblem(sunrealtype beta_val  = SUN_RCONST(10.0),
             sunrealtype gamma_val = SUN_RCONST(10.0))
    : beta(beta_val), gamma(gamma_val)
  {}

  // Get number of equations
  sunrealtype getBeta() const { return beta; }

  sunrealtype getGamma() const { return gamma; }

  int getNumEquations() const { return NEQ; }

  // Set initial conditions
  void setInitialConditions(N_Vector y)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    ydata[0]           = SUN_RCONST(1.0); // y1(0) = 1
    ydata[1]           = SUN_RCONST(1.0); // y2(0) = 1
    ydata[2]           = SUN_RCONST(1.0); // y3(0) = 1
  }

  // Right-hand side function (instance method)
  int computeRHS(sunrealtype t, N_Vector y, N_Vector ydot)
  {
    sunrealtype* ydata  = N_VGetArrayPointer(y);
    sunrealtype* dydata = N_VGetArrayPointer(ydot);

    sunrealtype y1 = ydata[0];
    sunrealtype y2 = ydata[1];
    sunrealtype y3 = ydata[2];

    // y1' = -10*y1 + beta*y2
    dydata[0] = -10.0 * y1 + beta * y2;

    // y2' = -beta*y1 - 10*y2
    dydata[1] = -beta * y1 - 10.0 * y2;

    // y3' = -gamma*y3
    dydata[2] = -gamma * y3;

    return 0; // Success
  }

  // Static wrapper for RHS (to pass to SUNDIALS)
  static int rhsWrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeRHS(t, y, ydot);
  }

  // Jacobian function (instance method)
  int computeJac(sunrealtype t, N_Vector y, SUNMatrix J)
  {
    // The Jacobian matrix is:
    // [-10    beta    0  ]
    // [-beta  -10     0  ]
    // [0      0      -gamma]

    // sunrealtype* Jdata = SUNDenseMatrix_Data(J);

    // Fill Jacobian (column-major order for SUNDIALS dense matrices)
    // Column 0
    SM_ELEMENT_D(J, 0, 0) = -10.0;
    SM_ELEMENT_D(J, 1, 0) = -beta;
    SM_ELEMENT_D(J, 2, 0) = 0.0;

    // Column 1
    SM_ELEMENT_D(J, 0, 1) = beta;
    SM_ELEMENT_D(J, 1, 1) = -10.0;
    SM_ELEMENT_D(J, 2, 1) = 0.0;

    // Column 2
    SM_ELEMENT_D(J, 0, 2) = 0.0;
    SM_ELEMENT_D(J, 1, 2) = 0.0;
    SM_ELEMENT_D(J, 2, 2) = -gamma;

    return 0; // Success
  }

  // Static wrapper for Jacobian (to pass to SUNDIALS)
  static int jacWrapper(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                        void* user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeJac(t, y, J);
  }

  // Return the eigenvalues of the Jacobian
  void computeEigenvalues(std::complex<sunrealtype>& lambda1,
                          std::complex<sunrealtype>& lambda2,
                          std::complex<sunrealtype>& lambda3) const
  {
    lambda1 = std::complex<sunrealtype>(-10, beta);   // -10 + beta i
    lambda2 = std::complex<sunrealtype>(-10, -beta);  // -10 - beta i
    lambda3 = std::complex<sunrealtype>(-gamma, 0.0); // -gamma
  }

  sunrealtype computeStiffnessRatio() const
  {
    std::complex<sunrealtype> lambda1, lambda2, lambda3;
    computeEigenvalues(lambda1, lambda2, lambda3);

    auto mag1 = std::abs(lambda1); // mag1 = mag2
    auto mag3 = std::abs(lambda3);

    auto min_mag = std::min(mag1, mag3);
    auto max_mag = std::max(mag1, mag3);

    if (min_mag > 0.0) { return max_mag / min_mag; }
    else
    {
      if (max_mag > 0.0)
      {
        return std::numeric_limits<sunrealtype>::infinity();
      }
      else
      {
        return 1.0;
      }
    }
  }
};

} // namespace shampine
} // namespace problems

#endif // SHAMPINE_HPP_
