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

#include <algorithm>
#include <array>
#include <complex>
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

#include <sundials/sundials_lapack_defs.h>

/* Interfaces to match 'sunrealtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgeev_f77 dgeev_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgeev_f77 sgeev_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif

namespace problems {
namespace rlc_ode {

struct rlc_eigenvalues
{
  std::array<std::complex<sunrealtype>, 4> eigenvalues{};
  std::array<std::array<std::complex<sunrealtype>, 4>, 4> rightEigenvectors{};
  std::array<std::array<std::complex<sunrealtype>, 4>, 4> leftEigenvectors{};
  sunrealtype min_mag;
  sunrealtype max_mag;
  int min_mag_idx;
  int max_mag_idx;

  void print(std::ostream& os = std::cout, bool printHeader = false, bool printLeftEigenvectors = false) const
  {
    auto print_complex = [&](const std::complex<sunrealtype>& z)
    {
      const auto re = static_cast<long double>(z.real());
      const auto im = static_cast<long double>(z.imag());
      os << re << (im < 0 ? " - " : " + ") << std::abs(im) << "i";
    };

    if (printHeader)
    {
      os << "eigenvalue, right_eigenvector";
      if (printLeftEigenvectors) { os << ", eigenvalue left_eigenvector"; }
      os << "\n";
    }
    os << std::scientific;
    os << std::setprecision(std::numeric_limits<sunrealtype>::digits10);

    for (int i = 0; i < 4; ++i)
    {
      print_complex(eigenvalues[i]);
      for (int k = 0; k < 4; ++k)
      {
        os << ", ";
        print_complex(rightEigenvectors[i][k]);
      }
      if (printLeftEigenvectors)
      {
        for (int k = 0; k < 4; ++k)
        {
          os << ", ";
          print_complex(leftEigenvectors[i][k]);
        }
      }
      os << "\n";
    }
  }
};

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

  void setVoltageSource(sunrealtype vs) { Vs = vs; }

  void printHeader(const std::vector<std::string>& labels, std::ostream& os = std::cout) const
  {
    os << std::setw(26) << "time";
    for (const auto& label : labels) {
      os << std::setw(26) << label;
    }
    os << "\n";
  }

  void printHeader(std::initializer_list<std::string> labels, std::ostream& os = std::cout) const
  {
    printHeader(std::vector<std::string>(labels), os);
  }

  void printState(sunrealtype t, N_Vector y, std::ostream& os = std::cout) const
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    os << std::scientific;
    os << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
    os << std::setw(26) << t;
    for (int i = 0; i < NEQ; ++i){ os << std::setw(26) << ydata[i]; }
    os << "\n";
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

    // Column 2 - d/dvC4
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

  // Compute eigenvalues of the Jacobian at a given state
  rlc_eigenvalues computeEigenvalues(bool computeLeftEigenvectors = false) const
  {
    constexpr int N = 4;

    // Convert to double column-major for LAPACK.
    sunrealtype MinvJ[16];
    // Column 0 - d/di2
    MinvJ[0] = (-(R1 + R2) - R3) / L2;
    MinvJ[1] = 1.0 / L2;
    MinvJ[2] = R3 / L2;
    MinvJ[3] = 0.0;
    // Column 1 - d/dvC3
    MinvJ[4] = -1.0 / C3;
    MinvJ[5] = 0.0;
    MinvJ[6] = 1.0 / C3;
    MinvJ[7] = 0.0;
    // Column 2 - d/di4
    MinvJ[8]  = R3 / L4;
    MinvJ[9]  = -1.0 / L4;
    MinvJ[10] = (-R3 - R4) / L4;
    MinvJ[11] = 1.0;
    // Column 2 - d/dvC4
    MinvJ[12] = 0.0;
    MinvJ[13] = 0.0;
    MinvJ[14] = -1.0 / C4;
    MinvJ[15] = 0.0;

    sunrealtype WR[N], WI[N];

    char JOBVL = computeLeftEigenvectors ? 'V' : 'N';
    char JOBVR = 'V';
    int LDA    = N;
    int LDVL   = N;
    int LDVR   = N;
    int INFO   = 0;

    // LAPACK outputs VL/VR in column-major, each is N x N when requested.
    std::vector<sunrealtype> VL_mat(computeLeftEigenvectors ? (N * N) : 1);
    std::vector<sunrealtype> VR_mat(N * N);

    // Workspace query
    int LWORK              = -1;
    sunrealtype WORK_QUERY = 0.0;
    xgeev_f77(&JOBVL, &JOBVR, (int*)&N, MinvJ, &LDA, WR, WI, VL_mat.data(),
              &LDVL, VR_mat.data(), &LDVR, &WORK_QUERY, &LWORK, &INFO);

    if (INFO != 0)
      throw std::runtime_error(
        "LAPACK workspace query failed (INFO=" + std::to_string(INFO) + ")");

    LWORK = std::max(1, static_cast<int>(WORK_QUERY));
    std::vector<sunrealtype> WORK(static_cast<size_t>(LWORK));

    // Compute eigenvalues/eigenvectors
    xgeev_f77(&JOBVL, &JOBVR, (int*)&N, MinvJ, &LDA, WR, WI, VL_mat.data(),
              &LDVL, VR_mat.data(), &LDVR, WORK.data(), &LWORK, &INFO);

    if (INFO != 0)
      throw std::runtime_error("GEEV failed (INFO=" + std::to_string(INFO) + ")");

    rlc_eigenvalues out{};

    // Eigenvalues
    for (int i = 0; i < N; ++i)
      out.eigenvalues[i] = {static_cast<sunrealtype>(WR[i]),
                            static_cast<sunrealtype>(WI[i])};

    auto get_col = [](const std::vector<sunrealtype>& M, int n, int col,
                      int row) -> sunrealtype
    {
      // M is column-major n x n
      return M[col * n + row];
    };

    // For a real eigenvalue WR[i] (WI[i]==0): the i-th right eigenvector is VR(:,i) (real).
    //
    // For a complex conjugate pair (WI[i]>0, WI[i+1]<0):
    //   right eigenvector for eigenvalue (WR[i] + i*WI[i]) is VR(:,i) + i*VR(:,i+1)
    //   right eigenvector for eigenvalue (WR[i] - i*WI[i]) is VR(:,i) - i*VR(:,i+1)

    // Right eigenvectors (complex reconstruction for conjugate pairs)
    for (int i = 0; i < N; ++i)
    {
      if (WI[i] == 0.0)
      {
        for (int r = 0; r < N; ++r)
          out.rightEigenvectors[i][r] = {static_cast<sunrealtype>(
                                           get_col(VR_mat, N, i, r)),
                                         static_cast<sunrealtype>(0.0)};
      }
      else if (WI[i] > 0.0)
      {
        // i and i+1 form a conjugate pair
        for (int r = 0; r < N; ++r)
        {
          const sunrealtype a             = get_col(VR_mat, N, i, r);
          const sunrealtype b             = get_col(VR_mat, N, i + 1, r);
          out.rightEigenvectors[i][r]     = {static_cast<sunrealtype>(a),
                                             static_cast<sunrealtype>(b)};
          out.rightEigenvectors[i + 1][r] = {static_cast<sunrealtype>(a),
                                             static_cast<sunrealtype>(-b)};
        }
      }
      // WI[i] < 0.0 is handled when we processed the prior WI[i-1] > 0.0
    }

    // Left eigenvectors (optional, same reconstruction rules applied to VL)
    if (computeLeftEigenvectors)
    {
      for (int i = 0; i < N; ++i)
      {
        if (WI[i] == 0.0)
        {
          for (int r = 0; r < N; ++r)
            out.leftEigenvectors[i][r] = {static_cast<sunrealtype>(
                                            get_col(VL_mat, N, i, r)),
                                          static_cast<sunrealtype>(0.0)};
        }
        else if (WI[i] > 0.0)
        {
          for (int r = 0; r < N; ++r)
          {
            const sunrealtype a            = get_col(VL_mat, N, i, r);
            const sunrealtype b            = get_col(VL_mat, N, i + 1, r);
            out.leftEigenvectors[i][r]     = {static_cast<sunrealtype>(a),
                                              static_cast<sunrealtype>(b)};
            out.leftEigenvectors[i + 1][r] = {static_cast<sunrealtype>(a),
                                              static_cast<sunrealtype>(-b)};
          }
        }
      }
    }

    auto result = std::minmax_element(out.eigenvalues.begin(),
                                      out.eigenvalues.end(),
                                      [](std::complex<sunrealtype> a,
                                         std::complex<sunrealtype> b)
                                      { return std::abs(a) < std::abs(b); });

    out.min_mag     = std::abs(*result.first);
    out.max_mag     = std::abs(*result.second);
    out.min_mag_idx = std::distance(out.eigenvalues.begin(), result.first);
    out.max_mag_idx = std::distance(out.eigenvalues.begin(), result.second);

    return out;
  }
};

} // namespace rlc_ode
} // namespace problems

#endif // RLC_ODE_HPP_
