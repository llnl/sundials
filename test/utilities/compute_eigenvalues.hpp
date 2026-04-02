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
 * ---------------------------------------------------------------------------*/

#ifndef COMPUTE_EIGENVALUES_HPP_
#define COMPUTE_EIGENVALUES_HPP_

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

// // 1) Form J = M^{-1} A without explicit inverse.
// ColMajorMatrix J = form_jacobian_from_mass(M, A);

// // 2) Run DGEEV on J with JOBVL='V', JOBVR='V' to get:
// //    wr, wi, VL, VR

// auto modes = build_ode_modes_from_dgeev(
//     n,
//     wr, wi,
//     VL.data(), n,
//     VR.data(), n);

// // 3) Pick multiple stiff modes.
// auto stiff_modes = select_stiff_modes_top_k(modes, 2);
// // or:
// // auto stiff_modes = select_stiff_modes_by_relative_decay(modes, 0.1);

// // 4) Score variables.
// auto scores_right = classify_stiff_vars_right_only_ode(modes, stiff_modes);
// auto scores_both  = classify_stiff_vars_left_right_ode(modes, stiff_modes);

// // 5) Rank variables: index 0 means y_1, index 1 means y_2, etc.
// auto rank_right = rank_variables_desc(scores_right);
// auto rank_both  = rank_variables_desc(scores_both);

struct eigenvalues
{
  std::array<std::complex<sunrealtype>, 4> eigenvalues{};
  std::array<std::array<std::complex<sunrealtype>, 4>, 4> rightEigenvectors{};
  std::array<std::array<std::complex<sunrealtype>, 4>, 4> leftEigenvectors{};
  int min_mag_idx;
  int max_mag_idx;
  sunrealtype min_mag;
  sunrealtype max_mag;
  sunrealtype stiffness_ratio;

  void print(std::ostream& os = std::cout, bool printHeader = false) const
  {
    // save original flags
    std::ios::fmtflags old_settings = std::cout.flags();

    auto print_complex = [&](const std::complex<sunrealtype>& z)
    {
      const auto re = static_cast<long double>(z.real());
      const auto im = static_cast<long double>(z.imag());
      os << re << (im < 0 ? " - " : " + ") << std::abs(im) << "i";
    };

    if (printHeader)
    {
      os << "eigenvalue, right eigenvector, left eigenvector\n";
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
      for (int k = 0; k < 4; ++k)
      {
        os << ", ";
        print_complex(leftEigenvectors[i][k]);
      }
      os << "\n";
    }

    // Restore original flags
    std::cout.flags(old_settings);
  }
};

// Compute eigenvalues of the Jacobian at a given state
eigenvalues computeEigenvalues(sunrealtype* MinvJ)
{
  constexpr int N = 4;

  sunrealtype WR[N], WI[N];

  char JOBVL = 'V';
  char JOBVR = 'V';
  int LDA    = N;
  int LDVL   = N;
  int LDVR   = N;
  int INFO   = 0;

  // LAPACK outputs VL/VR in column-major, each is N x N when requested.
  std::vector<sunrealtype> VL_mat(N * N);
  std::vector<sunrealtype> VR_mat(N * N);

  // Workspace query
  int LWORK              = -1;
  sunrealtype WORK_QUERY = 0.0;
  xgeev_f77(&JOBVL, &JOBVR, (int*)&N, MinvJ, &LDA, WR, WI, VL_mat.data(), &LDVL,
            VR_mat.data(), &LDVR, &WORK_QUERY, &LWORK, &INFO);

  if (INFO != 0)
    throw std::runtime_error(
      "LAPACK workspace query failed (INFO=" + std::to_string(INFO) + ")");

  LWORK = std::max(1, static_cast<int>(WORK_QUERY));
  std::vector<sunrealtype> WORK(static_cast<size_t>(LWORK));

  // Compute eigenvalues/eigenvectors
  xgeev_f77(&JOBVL, &JOBVR, (int*)&N, MinvJ, &LDA, WR, WI, VL_mat.data(), &LDVL,
            VR_mat.data(), &LDVR, WORK.data(), &LWORK, &INFO);

  if (INFO != 0)
    throw std::runtime_error("GEEV failed (INFO=" + std::to_string(INFO) + ")");

  eigenvalues out{};

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

  auto result = std::minmax_element(out.eigenvalues.begin(),
                                    out.eigenvalues.end(),
                                    [](std::complex<sunrealtype> a,
                                       std::complex<sunrealtype> b)
                                    { return std::abs(a) < std::abs(b); });

  out.min_mag         = std::abs(*result.first);
  out.max_mag         = std::abs(*result.second);
  out.min_mag_idx     = std::distance(out.eigenvalues.begin(), result.first);
  out.max_mag_idx     = std::distance(out.eigenvalues.begin(), result.second);
  out.stiffness_ratio = out.max_mag / out.min_mag;

  return out;
}

#endif // COMPUTE_EIGENVALUES_
