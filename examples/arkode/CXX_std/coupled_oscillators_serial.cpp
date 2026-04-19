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

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

/* --------------------------------------------------------------------------
 * Problem parameters (passed through user_data)
 * -------------------------------------------------------------------------- */

struct OscParams
{
  int N;    /* number of oscillators */
  double k; /* coupling strength     */
};

/* --------------------------------------------------------------------------
 * Index layout convention
 *
 * The state is a flat array of length 2*N:
 *   y[2*i + 0] = q_i   (position of oscillator i)
 *   y[2*i + 1] = p_i   (momentum of oscillator i)
 * -------------------------------------------------------------------------- */

static inline sunindextype Q(int i) { return 2 * i; }

static inline sunindextype P(int i) { return 2 * i + 1; }

/* --------------------------------------------------------------------------
 * RHS: coupled oscillators on a ring
 * -------------------------------------------------------------------------- */

int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto& params = *static_cast<OscParams*>(user_data);
  int N        = params.N;
  double k     = params.k;

  /* Flat pointer access */
  const sunrealtype* ys = N_VGetArrayPointer(y);
  sunrealtype* fs       = N_VGetArrayPointer(ydot);

  for (int i = 0; i < N; ++i)
  {
    int ip = (i + 1) % N;     /* periodic neighbor right */
    int im = (i - 1 + N) % N; /* periodic neighbor left  */

    /* dq_i/dt = p_i */
    fs[Q(i)] = ys[P(i)];

    /* dp_i/dt = -q_i + k*(q_{i+1} - 2*q_i + q_{i-1}) */
    fs[P(i)] = -ys[Q(i)] + k * (ys[Q(ip)] - 2.0 * ys[Q(i)] + ys[Q(im)]);
  }

  return 0;
}

/* --------------------------------------------------------------------------
 * Compute total energy: E = sum_i [ 0.5*p_i^2 + 0.5*q_i^2
 *                                   + 0.5*k*(q_{i+1} - q_i)^2 ]
 * -------------------------------------------------------------------------- */

sunrealtype totalEnergy(const sunrealtype* y, int N, double k)
{
  sunrealtype E = 0;
  for (int i = 0; i < N; ++i)
  {
    int ip = (i + 1) % N;
    E += 0.5 * y[P(i)] * y[P(i)];                               /* kinetic    */
    E += 0.5 * y[Q(i)] * y[Q(i)];                               /* on-site    */
    E += 0.5 * k * (y[Q(ip)] - y[Q(i)]) * (y[Q(ip)] - y[Q(i)]); /* coupling */
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

  const int neq = params.N * 2; /* total number of scalar equations */

  N_Vector y         = N_VNew_Serial(neq, sunctx);
  sunrealtype* ydata = N_VGetArrayPointer(y);

  /* Excite the first oscillator: q_0 = 1, all others zero. */
  N_VConst(0.0, y);
  ydata[Q(0)] = 1.0;

  void* ark_mem = ERKStepCreate(rhs, t0, y, sunctx);
  ARKodeSetUserData(ark_mem, &params);
  ARKodeSStolerances(ark_mem, 1.0e-6, 1.0e-14);

  sunrealtype E0 = totalEnergy(ydata, params.N, params.k);

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

    sunrealtype E = totalEnergy(ydata, params.N, params.k);
    std::printf("  %8.2f  %14.10f  %14.10f  %10.2e\n", tret, E, E0,
                std::abs((E - E0) / E0));

    tout += dt;
  }

  std::printf("\n  Final positions (q_i at t = %.1f):\n", tret);
  for (int i = 0; i < params.N; ++i)
    std::printf("    q_%d = %+.8f\n", i, ydata[Q(i)]);

  ARKodePrintAllStats(ark_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  ARKodeFree(&ark_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return 0;
}
