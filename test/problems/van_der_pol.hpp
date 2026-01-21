#ifndef VAN_DER_POL_HPP_
#define VAN_DER_POL_HPP_

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace problems {
namespace van_der_pol {

class ODEProblem
{
private:
  sunrealtype mu;           // Stiffness parameter
  static const int NEQ = 2; // Number of equations

public:
  // Constructor
  ODEProblem(sunrealtype mu_val = SUN_RCONST(10.0))
  : mu(mu_val) {}

  // Get parameters
  sunrealtype getMu() const { return mu; }
  int getNumEquations() const { return NEQ; }

  // Set initial conditions
  void setInitialConditions(N_Vector y)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    ydata[0]           = 2.0; // y1(0) = 2.0
    ydata[1]           = 0.0; // y2(0) = 0.0
  }

  // Right-hand side function (instance method)
  // Van der Pol equation: y1'' - mu*(1 - y1^2)*y1' + y1 = 0
  // Rewritten as first-order system:
  // y1' = y2
  // y2' = mu*(1 - y1^2)*y2 - y1
  int computeRHS(sunrealtype t, N_Vector y, N_Vector ydot)
  {
    sunrealtype* ydata  = N_VGetArrayPointer(y);
    sunrealtype* dydata = N_VGetArrayPointer(ydot);

    sunrealtype y1 = ydata[0];
    sunrealtype y2 = ydata[1];

    // y1' = y2
    dydata[0] = y2;

    // y2' = mu*(1 - y1^2)*y2 - y1
    dydata[1] = mu * (1.0 - y1 * y1) * y2 - y1;

    return 0; // Success
  }

  // Jacobian function (instance method)
  int computeJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    sunrealtype y1     = ydata[0];
    sunrealtype y2     = ydata[1];

    // The Jacobian matrix is:
    // [0                    1              ]
    // [-2*mu*y1*y2 - 1      mu*(1 - y1^2)  ]

    // Column 0
    SM_ELEMENT_D(J, 0, 0) = 0.0;
    SM_ELEMENT_D(J, 1, 0) = -2.0 * mu * y1 * y2 - 1.0;

    // Column 1
    SM_ELEMENT_D(J, 0, 1) = 1.0;
    SM_ELEMENT_D(J, 1, 1) = mu * (1.0 - y1 * y1);

    return 0; // Success
  }

  // Static wrapper for RHS (to pass to SUNDIALS)
  static int rhsWrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeRHS(t, y, ydot);
  }

  // Static wrapper for Jacobian (to pass to SUNDIALS)
  static int jacWrapper(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                        void* user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeJac(t, y, fy, J);
  }
};

} // namespace van_der_pol
} // namespace problems

#endif // VAN_DER_POL_HPP_
