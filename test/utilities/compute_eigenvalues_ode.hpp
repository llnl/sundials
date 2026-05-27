#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <sundials/sundials_lapack_defs.h>
#include <sunmatrix/sunmatrix_dense.h>

/* Interfaces to match 'sunrealtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgetrf_f77 dgetrf_f77
#define xgetrs_f77 dgetrs_f77
#define xgeev_f77  dgeev_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgetrf_f77 sgetrf_f77
#define xgetrs_f77 sgetrs_f77
#define xgeev_f77  sgeev_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif

struct OdeMode
{
  std::complex<sunrealtype> lambda{};
  std::vector<std::complex<sunrealtype>> vr; // right eigenvector
  std::vector<std::complex<sunrealtype>> vl; // left eigenvector
  int lapack_index  = -1;
  bool complex_pair = false;
};

struct VariableScores
{
  std::vector<sunrealtype> raw;
  std::vector<sunrealtype> normalized;
};

// =============================================================================

// -----------------------------------------------------------------------------
// 3) Select multiple stiff modes
// -----------------------------------------------------------------------------

enum class ModeWeightKind
{
  DecayRate,           // old behavior: max(0,-Re(lambda))^p
  EigenvalueMagnitude, // |lambda|^p
  IntrinsicRate,       // max(max(0,-Re(lambda)), |Im(lambda)|)^p
  StabilityFunction    // (1 / h_stab)^p
};

struct ModeWeightOptions
{
  ModeWeightKind kind = ModeWeightKind::DecayRate;

  // Exponent applied to the chosen base weight.
  sunrealtype power = 1.0;

  // If true, modes with Re(lambda) > stable_tol are ignored.
  bool ignore_unstable   = true;
  sunrealtype stable_tol = 1e-12;

  // Used only for kind == StabilityFunction.
  std::function<std::complex<sunrealtype>(std::complex<sunrealtype>)> stability_function;

  // Search controls for h_stab along the ray z = h * lambda.
  // If initial_h <= 0, a heuristic is used.
  sunrealtype initial_h = -1.0;
  sunrealtype max_h     = 1e12;
  int max_expand_steps  = 80;
  int max_bisect_steps  = 60;

  // Stability test: |R(z)| <= 1 + stability_tol
  sunrealtype stability_tol = 1e-10;

  // If h_stab is numerically zero, use min_h to avoid infinity.
  sunrealtype min_h = 1e-15;
};

static bool finite_complex(const std::complex<sunrealtype>& z)
{ return std::isfinite(z.real()) && std::isfinite(z.imag()); }

static sunrealtype stable_decay(const std::complex<sunrealtype>& lambda,
                                sunrealtype stable_tol = 1e-12)
{
  if (!finite_complex(lambda)) { return 0.0; }
  return (lambda.real() < -stable_tol) ? -lambda.real() : 0.0;
}

static sunrealtype intrinsic_rate(const std::complex<sunrealtype>& lambda,
                                  bool ignore_unstable,
                                  sunrealtype stable_tol = 1e-12)
{
  if (!finite_complex(lambda)) { return 0.0; }

  const sunrealtype decay =
    ignore_unstable ? ((lambda.real() < -stable_tol) ? -lambda.real() : 0.0)
                    : std::abs(lambda.real());

  const sunrealtype osc = std::abs(lambda.imag());

  return std::max(decay, osc);
}

// -----------------------------------------------------------------------------
// 4) Helpers for normalization
// -----------------------------------------------------------------------------

static sunrealtype l2_norm(const std::vector<std::complex<sunrealtype>>& x)
{
  sunrealtype s = 0.0;
  for (const auto& xi : x) { s += std::norm(xi); }
  return std::sqrt(s);
}

static void normalize_l2(std::vector<std::complex<sunrealtype>>& x)
{
  const sunrealtype nrm = l2_norm(x);
  if (nrm == 0.0) { throw std::runtime_error("Cannot normalize zero vector"); }
  for (auto& xi : x) { xi /= nrm; }
}

static std::complex<sunrealtype> dotc(
  const std::vector<std::complex<sunrealtype>>& x,
  const std::vector<std::complex<sunrealtype>>& y)
{
  if (x.size() != y.size())
  {
    throw std::runtime_error("Dimension mismatch in dotc");
  }

  std::complex<sunrealtype> s = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) { s += std::conj(x[i]) * y[i]; }
  return s;
}

static VariableScores finalize_scores(const std::vector<sunrealtype>& raw)
{
  VariableScores out;
  out.raw        = raw;
  out.normalized = raw;

  const sunrealtype sum = std::accumulate(raw.begin(), raw.end(), 0.0);
  if (sum > 0.0)
  {
    for (sunrealtype& s : out.normalized) { s /= sum; }
  }
  return out;
}

// Estimate the largest stable step along z = h * lambda such that
// |R(h*lambda)| <= 1 + tol.
//
// Assumes the stable set along this ray is the interval [0, h_stab],
// which is the usual case for standard one-step methods near the origin.
static sunrealtype estimate_stable_step_on_ray(
  const std::complex<sunrealtype>& lambda, const ModeWeightOptions& opts)
{
  if (!opts.stability_function)
  {
    throw std::runtime_error("stability_function is not set");
  }

  if (!finite_complex(lambda)) { return 0.0; }

  if (std::abs(lambda) == 0.0)
  {
    return std::numeric_limits<sunrealtype>::infinity();
  }

  if (opts.ignore_unstable && lambda.real() > opts.stable_tol) { return 0.0; }

  auto stable_at = [&](sunrealtype h) -> bool
  {
    const std::complex<sunrealtype> z  = h * lambda;
    const std::complex<sunrealtype> Rz = opts.stability_function(z);
    const sunrealtype mag              = std::abs(Rz);
    return std::isfinite(mag) && mag <= 1.0 + opts.stability_tol;
  };

  if (!stable_at(0.0)) { return 0.0; }

  sunrealtype h_lo = 0.0;
  sunrealtype h_hi = opts.initial_h;

  if (h_hi <= 0.0) { h_hi = 1.0 / std::max(1.0, std::abs(lambda)); }
  if (h_hi <= 0.0) { h_hi = 1.0; }

  int expand_count = 0;
  while (expand_count < opts.max_expand_steps && h_hi < opts.max_h &&
         stable_at(h_hi))
  {
    h_lo = h_hi;
    h_hi = std::min(2.0 * h_hi, opts.max_h);
    ++expand_count;
  }

  if (stable_at(h_hi)) { return std::numeric_limits<sunrealtype>::infinity(); }

  for (int it = 0; it < opts.max_bisect_steps; ++it)
  {
    const sunrealtype h_mid = 0.5 * (h_lo + h_hi);
    if (stable_at(h_mid)) { h_lo = h_mid; }
    else
    {
      h_hi = h_mid;
    }
  }

  return h_lo;
}

static sunrealtype compute_mode_weight(const std::complex<double>& lambda,
                                       const ModeWeightOptions& opts)
{
  if (!finite_complex(lambda)) { return 0.0; }

  if (opts.ignore_unstable && lambda.real() > opts.stable_tol) { return 0.0; }

  switch (opts.kind)
  {
  case ModeWeightKind::DecayRate:
  {
    const double d = stable_decay(lambda, opts.stable_tol);
    return (d > 0.0) ? std::pow(d, opts.power) : 0.0;
  }

  case ModeWeightKind::EigenvalueMagnitude:
  {
    const double a = std::abs(lambda);
    return (a > 0.0) ? std::pow(a, opts.power) : 0.0;
  }

  case ModeWeightKind::IntrinsicRate:
  {
    const double r = intrinsic_rate(lambda, opts.ignore_unstable,
                                    opts.stable_tol);
    return (r > 0.0) ? std::pow(r, opts.power) : 0.0;
  }

  case ModeWeightKind::StabilityFunction:
  {
    const double h_stab = estimate_stable_step_on_ray(lambda, opts);

    // Infinite stable step => this mode does not impose a stability limit.
    if (!std::isfinite(h_stab)) { return 0.0; }

    // Very small or zero stable step => very restrictive mode.
    const double h_eff = std::max(h_stab, opts.min_h);
    return std::pow(1.0 / h_eff, opts.power);
  }
  }

  return 0.0;
}

// Optional helper: choose the k modes with the largest weight.
// This lets you select multiple "stiff" modes using the same criterion
// that will later be used in scoring.
std::vector<int> select_modes_top_k_by_weight(const std::vector<OdeMode>& modes,
                                              int k,
                                              const ModeWeightOptions& opts)
{
  std::vector<std::pair<double, int>> scored;

  for (int i = 0; i < static_cast<int>(modes.size()); ++i)
  {
    const double w = compute_mode_weight(modes[i].lambda, opts);
    if (w > 0.0) { scored.emplace_back(w, i); }
  }

  std::sort(scored.begin(), scored.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  if (k > 0 && static_cast<int>(scored.size()) > k) { scored.resize(k); }

  std::vector<int> out;
  out.reserve(scored.size());
  for (const auto& p : scored) { out.push_back(p.second); }
  return out;
}

// =============================================================================

void printModes(const std::vector<OdeMode>& modes, std::ostream& os = std::cout,
                bool printHeader = false)
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

  for (const auto& m : modes)
  {
    print_complex(m.lambda);
    for (const auto& r : m.vr)
    {
      os << ", ";
      print_complex(r);
    }
    for (const auto& l : m.vl)
    {
      os << ", ";
      print_complex(l);
    }
    os << "\n";
    if (m.complex_pair)
    {
      print_complex(std::conj(m.lambda));
      for (const auto& r : m.vr)
      {
        os << ", ";
        print_complex(std::conj(r));
      }
      for (const auto& l : m.vl)
      {
        os << ", ";
        print_complex(std::conj(l));
      }
      os << "\n";
    }
  }

  // Restore original flags
  std::cout.flags(old_settings);
}

// -----------------------------------------------------------------------------
// 1) Form J = M^{-1} A without explicitly forming M^{-1}
// -----------------------------------------------------------------------------
SUNMatrix form_jacobian_from_mass(SUNMatrix M, SUNMatrix A)
{
  int n          = SUNDenseMatrix_Rows(M);
  SUNMatrix M_lu = SUNMatClone(M);
  SUNMatrix J    = SUNMatClone(A); // overwritten with solution of M * J = A
  std::vector<int> ipiv(n);
  int info = 0;

  SUNMatCopy(M, M_lu);
  SUNMatCopy(A, J);

  xgetrf_f77(&n, &n, SUNDenseMatrix_Data(M_lu), &n, ipiv.data(), &info);
  if (info < 0) { throw std::runtime_error("DGETRF: illegal argument"); }
  if (info > 0) { throw std::runtime_error("DGETRF: M is singular"); }

  const char trans = 'N';
  const int nrhs   = n; // solve for all columns of A at once
  xgetrs_f77(&trans, &n, &nrhs, SUNDenseMatrix_Data(M_lu), &n, ipiv.data(),
             SUNDenseMatrix_Data(J), &n, &info);

  if (info != 0) { throw std::runtime_error("DGETRS failed"); }

  return J;
}

// -----------------------------------------------------------------------------
// 2) Unpack DGEEV output into unique complex modes
//
// DGEEV returns:
//   real eigenvalue:   lambda_j = WR[j]
//   complex pair:      lambda_j = WR[j] + i*WI[j]
//                      lambda_{j+1} = WR[j] - i*WI[j]
//
// and the eigenvectors for a complex pair are stored as:
//   v_j     = VR(:,j) + i*VR(:,j+1)
//   v_{j+1} = VR(:,j) - i*VR(:,j+1)
//
// Same for VL.
// We keep only the member with positive imaginary part to avoid double counting.
// -----------------------------------------------------------------------------
static std::vector<std::complex<sunrealtype>> unpack_dgeev_vector(
  const std::vector<sunrealtype>& V, int ldv, int n, int j, bool complex_pair)
{
  std::vector<std::complex<sunrealtype>> out(n);
  if (!complex_pair)
  {
    for (int i = 0; i < n; ++i)
    {
      out[i] = std::complex<sunrealtype>(V[i + j * ldv], 0.0);
    }
  }
  else
  {
    for (int i = 0; i < n; ++i)
    {
      out[i] = std::complex<sunrealtype>(V[i + j * ldv], V[i + (j + 1) * ldv]);
    }
  }
  return out;
}

std::vector<OdeMode> build_ode_modes_from_dgeev(SUNMatrix A,
                                                sunrealtype imag_tol = 1e-14)
{
  const int n        = SUNDenseMatrix_Rows(A);
  sunrealtype* Adata = SUNDenseMatrix_Data(A);

  std::vector<sunrealtype> wr(n);
  std::vector<sunrealtype> wi(n);

  char JOBVL = 'V';
  char JOBVR = 'V';
  int LDA    = n;
  int ldvl   = n;
  int ldvr   = n;
  int INFO   = 0;

  // LAPACK outputs VL/VR in column-major, each is N x N when requested.
  std::vector<sunrealtype> vl(n * n);
  std::vector<sunrealtype> vr(n * n);

  // Workspace query
  int LWORK              = -1;
  sunrealtype WORK_QUERY = 0.0;
  xgeev_f77(&JOBVL, &JOBVR, (int*)&n, Adata, &LDA, wr.data(), wi.data(),
            vl.data(), &ldvl, vr.data(), &ldvr, &WORK_QUERY, &LWORK, &INFO);

  if (INFO != 0)
    throw std::runtime_error(
      "LAPACK workspace query failed (INFO=" + std::to_string(INFO) + ")");

  LWORK = std::max(1, static_cast<int>(WORK_QUERY));
  std::vector<sunrealtype> WORK(static_cast<size_t>(LWORK));

  // Compute eigenvalues/eigenvectors
  xgeev_f77(&JOBVL, &JOBVR, (int*)&n, Adata, &LDA, wr.data(), wi.data(),
            vl.data(), &ldvl, vr.data(), &ldvr, WORK.data(), &LWORK, &INFO);

  if (INFO != 0)
    throw std::runtime_error("GEEV failed (INFO=" + std::to_string(INFO) + ")");

  std::vector<OdeMode> modes;
  modes.reserve(n);

  for (int j = 0; j < n; ++j)
  {
    if (wi[j] < -imag_tol)
    {
      continue; // second half of conjugate pair
    }

    const bool is_real = std::abs(wi[j]) <= imag_tol;
    const bool is_pair = !is_real;

    if (is_pair && j + 1 >= n)
    {
      throw std::runtime_error("Malformed DGEEV output");
    }

    OdeMode m;
    m.lambda       = std::complex<sunrealtype>(wr[j], wi[j]);
    m.vr           = unpack_dgeev_vector(vr, ldvr, n, j, is_pair);
    m.vl           = unpack_dgeev_vector(vl, ldvl, n, j, is_pair);
    m.lapack_index = j;
    m.complex_pair = is_pair;

    modes.push_back(std::move(m));
  }

  return modes;
}

// -----------------------------------------------------------------------------
// 5) Classification approach 1: right eigenvectors only
//
// Per selected mode i:
//   normalize v_i so ||v_i||_2 = 1
//   per-variable score = |v_i(j)|^2
//
// Aggregate over multiple stiff modes:
//   score_j = sum_i weight_i * |v_i(j)|^2
// -----------------------------------------------------------------------------
VariableScores classify_stiff_vars_right_only_ode(
  const std::vector<OdeMode>& modes, const std::vector<int>& selected_modes,
  const ModeWeightOptions& weight_opts)
{
  int n = -1;
  for (const auto& m : modes)
  {
    if (!m.vr.empty())
    {
      n = static_cast<int>(m.vr.size());
      break;
    }
  }
  if (n < 0) { throw std::runtime_error("No right eigenvectors available"); }

  std::vector<sunrealtype> raw(n, 0.0);

  for (int idx : selected_modes)
  {
    if (idx < 0 || idx >= static_cast<int>(modes.size()))
    {
      throw std::runtime_error("selected mode index out of range");
    }

    const auto& mode = modes[idx];
    if (mode.vr.empty())
    {
      throw std::runtime_error("Missing right eigenvector");
    }

    const double weight = compute_mode_weight(mode.lambda, weight_opts);
    if (weight == 0.0) { continue; }

    std::vector<std::complex<double>> v = mode.vr;
    normalize_l2(v);

    for (int j = 0; j < n; ++j) { raw[j] += weight * std::norm(v[j]); }
  }

  return finalize_scores(raw);
}

// -----------------------------------------------------------------------------
// 6) Classification approach 2: left + right eigenvectors
//
// Standard ODE case:
//   J v_i = lambda_i v_i
//   w_i^H J = lambda_i w_i^H
//
// Normalize each mode so that:
//   w_i^H v_i = 1
//
// Then use a participation-like per-variable profile:
//   p_j = |v_i(j) * w_i(j)|
//
// Normalize p over j for each mode, then aggregate over selected modes.
// -----------------------------------------------------------------------------
VariableScores classify_stiff_vars_left_right_ode(
  const std::vector<OdeMode>& modes, const std::vector<int>& selected_modes,
  const ModeWeightOptions& weight_opts, double biorth_tol = 1e-12)
{
  int n = -1;
  for (const auto& m : modes)
  {
    if (!m.vr.empty())
    {
      n = static_cast<int>(m.vr.size());
      break;
    }
  }
  if (n < 0) { throw std::runtime_error("No eigenvectors available"); }

  std::vector<sunrealtype> raw(n, 0.0);

  for (int idx : selected_modes)
  {
    if (idx < 0 || idx >= static_cast<int>(modes.size()))
    {
      throw std::runtime_error("stiff mode index out of range");
    }

    const auto& mode = modes[idx];
    if (mode.vr.empty() || mode.vl.empty())
    {
      throw std::runtime_error("Both left and right eigenvectors are required");
    }

    const double weight = compute_mode_weight(mode.lambda, weight_opts);
    if (weight == 0.0) { continue; }

    std::vector<std::complex<sunrealtype>> v = mode.vr;
    std::vector<std::complex<sunrealtype>> w = mode.vl;

    normalize_l2(v);

    const std::complex<sunrealtype> gamma = dotc(w, v); // w^H v
    if (std::abs(gamma) < biorth_tol)
    {
      throw std::runtime_error("w^H v is too small for stable normalization");
    }

    const std::complex<sunrealtype> scale_w = 1.0 / std::conj(gamma);
    for (auto& wi : w) { wi *= scale_w; }

    std::vector<sunrealtype> p(n, 0.0);
    sunrealtype psum = 0.0;
    for (int j = 0; j < n; ++j)
    {
      p[j] = std::abs(v[j] * w[j]);
      psum += p[j];
    }
    if (psum <= 0.0) { continue; }

    for (int j = 0; j < n; ++j) { raw[j] += weight * (p[j] / psum); }
  }

  return finalize_scores(raw);
}

std::vector<int> rank_variables_desc(const VariableScores& scores)
{
  const std::vector<sunrealtype>& s =
    !scores.normalized.empty() ? scores.normalized : scores.raw;

  std::vector<int> idx(s.size());
  std::iota(idx.begin(), idx.end(), 0);

  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return s[a] > s[b]; });

  return idx;
}
