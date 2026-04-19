/* -----------------------------------------------------------------------------
 * Solves a system of N coupled oscillators arranged on a ring:
 *
 *   dq_i/dt =  p_i
 *   dp_i/dt = -q_i + k*(q_{i+1} - 2*q_i + q_{i-1})
 *
 * with periodic boundary conditions (indices mod N) and coupling k.
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <mdspan>
#include <vector>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_stdmdspan.hpp>

/* --------------------------------------------------------------------------
 * Type aliases
 *
 * The state is an N x 2 matrix (N oscillators, 2 components each).
 * layout_right gives row-major storage: state[i, 0] and state[i, 1] are
 * contiguous for each oscillator.
 * -------------------------------------------------------------------------- */

using Span2D    = std::mdspan<sunrealtype, std::dextents<size_t, 2>>;
using MdSpanVec = sundials::example::MdSpanVector<Span2D>;

/* --------------------------------------------------------------------------
 * Problem parameters (passed through user_data)
 * -------------------------------------------------------------------------- */

struct OscParams
{
  int N;    /* number of oscillators */
  double k; /* coupling strength     */
};

/* --------------------------------------------------------------------------
 * RHS: coupled oscillators on a ring
 * -------------------------------------------------------------------------- */

int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto& params = *static_cast<OscParams*>(user_data);
  int N        = params.N;
  double k     = params.k;

  /* Recover the 2D mdspan from the N_Vector. */
  const auto& ys = MdSpanVec::extractConst(y)->span(); /* N x 2, read-only  */
  auto& fs       = MdSpanVec::extract(ydot)->span();   /* N x 2, write      */

  for (int i = 0; i < N; ++i)
  {
    int ip = (i + 1) % N;     /* periodic neighbor right */
    int im = (i - 1 + N) % N; /* periodic neighbor left  */

    /* dq_i/dt = p_i */
    fs[i, 0] = ys[i, 1];

    /* dp_i/dt = -q_i + k*(q_{i+1} - 2*q_i + q_{i-1}) */
    fs[i, 1] = -ys[i, 0] + k * (ys[ip, 0] - 2.0 * ys[i, 0] + ys[im, 0]);
  }

  return 0;
}

/* --------------------------------------------------------------------------
 * Compute total energy: E = sum_i [ 0.5*p_i^2 + 0.5*q_i^2
 *                                   + 0.5*k*(q_{i+1} - q_i)^2 ]
 * -------------------------------------------------------------------------- */

sunrealtype totalEnergy(const Span2D& s, int N, double k)
{
  sunrealtype E = 0;
  for (int i = 0; i < N; ++i)
  {
    int ip = (i + 1) % N;
    E += 0.5 * s[i, 1] * s[i, 1];                               /* kinetic    */
    E += 0.5 * s[i, 0] * s[i, 0];                               /* on-site    */
    E += 0.5 * k * (s[ip, 0] - s[i, 0]) * (s[ip, 0] - s[i, 0]); /* coupling */
  }
  return E;
}

/* --------------------------------------------------------------------------
 * Main
 * -------------------------------------------------------------------------- */

int main()
{
  /* Parameters */
  OscParams params;
  params.N = 8;    /* 8 oscillators on a ring */
  params.k = 0.25; /* moderate coupling       */

  const sunrealtype t0 = 0.0;
  const sunrealtype tf = 50.0;
  const sunrealtype dt = 5.0;

  SUNContext sunctx = nullptr;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Initial condition: N x 2 matrix via mdspan
   *
   * Storage holds the raw data; the mdspan provides 2D (i, component) indexing
   * over it; the MdSpanVector wraps the mdspan as an N_Vector. All three refer
   * to the same memory.
   */

  const int neq = params.N * 2; /* total number of scalar equations */

  std::vector<sunrealtype> storage(neq, 0.0);
  Span2D state(storage.data(), params.N, 2);

  /* Excite the first oscillator: q_0 = 1, all others zero. */
  state[0, 0] = 1.0;

  /* Wrap the mdspan as a SUNDIALS vector (view - no copy). */
  MdSpanVec y(state, sunctx);

  void* ark_mem = ERKStepCreate(rhs, t0, y, sunctx);
  ARKodeSetUserData(ark_mem, &params);
  ARKodeSStolerances(ark_mem, 1.0e-6, 1.0e-14);

  sunrealtype E0 = totalEnergy(state, params.N, params.k);

  std::printf("  %8s  %14s  %14s  %10s\n", "t", "energy", "E0", "|dE/E0|");
  std::printf("  --------------------------------------------------------\n");
  std::printf("  %8.2f  %14.10f  %14.10f  %10.2e\n", t0, E0, E0, 0.0);

  sunrealtype tret = t0;
  sunrealtype tout = t0 + dt;

  while (tout <= tf + 0.5 * dt)
  {
    int err = ARKodeEvolve(ark_mem, tout, y, &tret, ARK_NORMAL);
    if (err < 0)
    {
      std::fprintf(stderr, "ARKodeEvolve failed (err = %d)\n", err);
      break;
    }

    sunrealtype E = totalEnergy(state, params.N, params.k);
    std::printf("  %8.2f  %14.10f  %14.10f  %10.2e\n", tret, E, E0,
                std::abs((E - E0) / E0));

    tout += dt;
  }

  std::printf("\n  Final positions (q_i at t = %.1f):\n", tret);
  for (int i = 0; i < params.N; ++i)
    std::printf("    q_%d = %+.8f\n", i, state[i, 0]);

  ARKodePrintAllStats(ark_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  ARKodeFree(&ark_mem);
  SUNContext_Free(&sunctx);

  return 0;
}
