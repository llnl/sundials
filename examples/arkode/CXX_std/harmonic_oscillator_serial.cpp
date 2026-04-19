/* -----------------------------------------------------------------------------
 * Solves the 2D harmonic oscillator:
 *
 *   dy1/dt =  y2        y1(0) = 1
 *   dy2/dt = -y1        y2(0) = 0
 *
 * Exact solution: y1(t) = cos(t),  y2(t) = sin(t).
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

/* --------------------------------------------------------------------------
 * Right-hand-side function
 * -------------------------------------------------------------------------- */

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* /* user_data */)
{
  auto yd = N_VGetArrayPointer(y);
  auto fd = N_VGetArrayPointer(ydot);

  /* dy1/dt =  y2
   * dy2/dt = -y1 */
  fd[0] = yd[1];
  fd[1] = -yd[0];

  return 0;
}

/* --------------------------------------------------------------------------
 * Main
 * -------------------------------------------------------------------------- */

int main()
{
  SUNContext sunctx;
  int err = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (err != 0)
  {
    std::fprintf(stderr, "SUNContext_Create failed\n");
    return 1;
  }

  const sunindextype neq = 2;    /* number of equations */
  const sunrealtype t0   = 0.0;  /* initial time        */
  const sunrealtype tf   = 10.0; /* final time          */
  const sunrealtype dt   = 1.0;  /* output interval     */

  N_Vector y = N_VNew_Serial(neq, sunctx);
  auto ydata = N_VGetArrayPointer(y);
  ydata[0]   = 1.0; /* y1(0) = cos(0) = 1 */
  ydata[1]   = 0.0; /* y2(0) = sin(0) = 0 */

  void* arkode_mem = ERKStepCreate(rhs, t0, y, sunctx);
  if (arkode_mem == nullptr)
  {
    std::fprintf(stderr, "ERKStepCreate failed\n");
    SUNContext_Free(&sunctx);
    return 1;
  }

  sunrealtype reltol = 1.0e-8;
  sunrealtype abstol = 1.0e-12;
  ARKodeSStolerances(arkode_mem, reltol, abstol);

  std::printf("  %8s  %14s  %14s  %14s  %14s  %10s\n", "t", "y1", "y2",
              "y1 exact", "y2 exact", "error");
  std::printf("  -------------------------------------------"
              "---------------------------------------------------\n");

  /* Print initial condition */
  std::printf("  %8.4f  %14.10f  %14.10f  %14.10f  %14.10f  %10.2e\n", t0,
              ydata[0], ydata[1], std::cos(t0), -std::sin(t0),
              std::hypot(ydata[0] - std::cos(t0), ydata[1] + std::sin(t0)));

  sunrealtype tret = t0;
  sunrealtype tout = t0 + dt;

  while (tout <= tf + 0.5 * dt)
  {
    err = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
    if (err < 0)
    {
      std::fprintf(stderr, "ARKodeEvolve failed at t = %g (err = %d)\n", tout,
                   err);
      break;
    }

    /* Compare against exact solution */
    sunrealtype y1_exact = std::cos(tret);
    sunrealtype y2_exact = -std::sin(tret);
    sunrealtype error    = std::hypot(ydata[0] - y1_exact, ydata[1] - y2_exact);

    std::printf("  %8.4f  %14.10f  %14.10f  %14.10f  %14.10f  %10.2e\n", tret,
                ydata[0], ydata[1], y1_exact, y2_exact, error);

    tout += dt;
  }

  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  ARKodeFree(&arkode_mem);
  SUNContext_Free(&sunctx);

  return 0;
}
