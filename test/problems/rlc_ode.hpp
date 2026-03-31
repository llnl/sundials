/* -----------------------------------------------------------------------------
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
 * SPDX-License-Identifier: BSD-3-Clause SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This header defines and ODE model for an RLC circuit (resistor (R), inductor
 * (L), and capacitor (C).
 *
 * L2 *  i2' = Vs - (R1 + R2) * i2 - R3 * (i2 - i4) - vC3
 * C3 * vC3' = i2 - i4
 * L4 *  i4' = R3 * (i2 - i4) + vC3 - R4 * i4 - vC4
 * C4 * vC4' = i4
 *
 * Differential states:
 *  i2  = current into the V2 node
 *  vC3 = voltage across capacitor C3
 *  i4  = current out of the V2 node
 *  vC4 = voltage across capacitor C4
 *
 * Independent input:
 *  Vs = voltage source
 *
 * Problem constants:
 *  R1 = resistance for resistor R1
 *  R2 = resistance for resistor R2
 *  R3 = resistance for resistor R3
 *  R4 = resistance for resistor R4
 *  L2 = inductance for inductor L2
 *  L4 = inductance for inductor L4
 *  C3 = capacitance for capacitor C3
 *  C4 = capacitance for capacitor C4
 * ---------------------------------------------------------------------------*/

#ifndef RLC_ODE_HPP_
#define RLC_ODE_HPP_

#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace problems {
namespace rlc_ode {

class ODEProblem
{
private:
  sunrealtype R1, R2, R3, R4; // Resistance
  sunrealtype L2, L4;         // Inductance
  sunrealtype C3, C4;         // Capacitance
  sunrealtype Vs;             // Voltage source
  static const int NEQ = 4;   // Number of equations

public:
  // Constructor
  ODEProblem(sunrealtype r1 = 10.0, sunrealtype r2 = 10.0, sunrealtype r3 = 2.0,
             sunrealtype r4 = 20.0, sunrealtype l2 = 10.0, sunrealtype l4 = 10.0,
             sunrealtype c3 = 10.0, sunrealtype c4 = 50.0, sunrealtype vs = 12)
    : R1(r1), R2(r2), R3(r3), R4(r4), L2(l2), L4(l4), C3(c3), C4(c4), Vs(vs)
  {}

  int getNumEquations() const { return NEQ; }

  sunrealtype getR1() { return R1; }

  sunrealtype getR2() { return R2; }

  sunrealtype getR3() { return R3; }

  sunrealtype getR4() { return R4; }

  sunrealtype getL2() { return L2; }

  sunrealtype getL4() { return L4; }

  sunrealtype getC3() { return C3; }

  sunrealtype getC4() { return C4; }

  sunrealtype getVs() { return Vs; }

  void setVoltageSource(sunrealtype vs) { Vs = vs; }

  void printHeader(const std::vector<std::string>& labels,
                   std::ostream& os = std::cout) const
  {
    os << std::setw(26) << "time";
    for (const auto& label : labels) { os << std::setw(26) << label; }
    os << "\n";
  }

  void printHeader(std::initializer_list<std::string> labels,
                   std::ostream& os = std::cout) const
  { printHeader(std::vector<std::string>(labels), os); }

  void printState(sunrealtype t, N_Vector y, std::ostream& os = std::cout) const
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);

    std::ios::fmtflags old_settings = std::cout.flags(); // save original flags
    os << std::scientific;
    os << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
    os << std::setw(26) << t;
    for (int i = 0; i < NEQ; ++i) { os << std::setw(26) << ydata[i]; }
    os << "\n";
    std::cout.flags(old_settings); // Restore original flags
  }

  // Mass matrix function (instance method)
  int computeMassDense(SUNMatrix M)
  {
    SUNMatZero(M);
    SM_ELEMENT_D(M, 0, 0) = L2;
    SM_ELEMENT_D(M, 1, 1) = C3;
    SM_ELEMENT_D(M, 2, 2) = L4;
    SM_ELEMENT_D(M, 3, 3) = C4;

    return 0; // Success
  }

  // Static wrapper for Mass matrix function (to pass to SUNDIALS)
  static int massDenseWrapper(sunrealtype t, SUNMatrix M, void* user_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeMassDense(M);
  }

  // Right-hand side function (instance method)
  int computeRHS(sunrealtype t, N_Vector y, N_Vector ydot)
  {
    sunrealtype* ydata  = N_VGetArrayPointer(y);
    sunrealtype* dydata = N_VGetArrayPointer(ydot);

    sunrealtype i2  = ydata[0];
    sunrealtype vC3 = ydata[1];
    sunrealtype i4  = ydata[2];
    sunrealtype vC4 = ydata[3];

    dydata[0] = Vs - (R1 + R2) * i2 - R3 * (i2 - i4) - vC3;
    dydata[1] = i2 - i4;
    dydata[2] = R3 * (i2 - i4) + vC3 - R4 * i4 - vC4;
    dydata[3] = i4;

    return 0; // Success
  }

  // Static wrapper for RHS (to pass to SUNDIALS)
  static int rhsWrapper(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeRHS(t, y, ydot);
  }

  // Jacobian function (instance method)
  int computeJacDense(SUNMatrix J)
  {
    // Column 0 - d/di2
    SM_ELEMENT_D(J, 0, 0) = -(R1 + R2) - R3;
    SM_ELEMENT_D(J, 1, 0) = 1.0;
    SM_ELEMENT_D(J, 2, 0) = R3;
    SM_ELEMENT_D(J, 3, 0) = 0.0;

    // Column 1 - d/dvC3
    SM_ELEMENT_D(J, 0, 1) = -1.0;
    SM_ELEMENT_D(J, 1, 1) = 0.0;
    SM_ELEMENT_D(J, 2, 1) = 1.0;
    SM_ELEMENT_D(J, 3, 1) = 0.0;

    // Column 2 - d/di4
    SM_ELEMENT_D(J, 0, 2) = R3;
    SM_ELEMENT_D(J, 1, 2) = -1.0;
    SM_ELEMENT_D(J, 2, 2) = -R3 - R4;
    SM_ELEMENT_D(J, 3, 2) = 1.0;

    // Column 3 - d/dvC4
    SM_ELEMENT_D(J, 0, 3) = 0.0;
    SM_ELEMENT_D(J, 1, 3) = 0.0;
    SM_ELEMENT_D(J, 2, 3) = -1.0;
    SM_ELEMENT_D(J, 3, 3) = 0.0;

    return 0; // Success
  }

  // Static wrapper for Jacobian (to pass to SUNDIALS)
  static int jacWrapperDense(sunrealtype t, N_Vector y, N_Vector fy,
                             SUNMatrix J, void* user_data, N_Vector tmp1,
                             N_Vector tmp2, N_Vector tmp3)
  {
    ODEProblem* problem = static_cast<ODEProblem*>(user_data);
    return problem->computeJacDense(J);
  }
};

} // namespace rlc_ode
} // namespace problems

#endif // RLC_ODE_HPP_
