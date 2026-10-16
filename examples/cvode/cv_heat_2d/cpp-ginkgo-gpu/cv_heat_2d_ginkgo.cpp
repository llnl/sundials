/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
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
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a simple anisotropic 2D heat equation,
 *
 *   u_t = kx u_xx + ky u_yy + b,
 *
 * for t in [0, 1] and (x,y) in [0, 1]^2, with initial conditions
 *
 *   u(0,x,y) = sin^2(pi x) sin^2(pi y) + 1,
 *
 * stationary boundary conditions
 *
 *   u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 *
 * and the heat source
 *
 *   b(t,x,y) = -2 pi sin^2(pi x) sin^2(pi y) sin(pi t) cos(pi t)
 *              - kx 2 pi^2 (cos^2(pi x) - sin^2(pi x)) sin^2(pi y) cos^2(pi t)
 *              - ky 2 pi^2 (cos^2(pi y) - sin^2(pi y)) sin^2(pi x) cos^2(pi t).
 *
 * Under this setup, the problem has the analytical solution
 *
 *    u(t,x,y) = sin^2(pi x) sin^2(pi y) cos^2(pi t) + 1.
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid. The
 * problem is advanced in time with BDF methods using an inexact Newton method
 * paired with the CG linear solver from Ginkgo. Several command line options
 * are available to change the problem parameters and CVODE settings. Use the
 * flag
 * --help for more information.
 * ---------------------------------------------------------------------------*/

// Include user data structure and utility functions for this problem
#include "cv_heat_2d_ginkgo.hpp"

#include <sunmatrix/sunmatrix_ginkgo.hpp>

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#define N_VNew_GPU                   N_VNew_Cuda
#define N_VGetDeviceArrayPointer_GPU N_VGetDeviceArrayPointer_Cuda
#define N_VCopyFromDevice_GPU        N_VCopyFromDevice_Cuda
#define GPU_DEVICE_SYNCHRONIZE       cudaDeviceSynchronize
#define GPU_GET_LAST_ERROR           cudaGetLastError
#define GPU_SUCCESS                  cudaSuccess
#define GPU_GET_ERROR_NAME           cudaGetErrorName
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#define N_VNew_GPU                   N_VNew_Hip
#define N_VGetDeviceArrayPointer_GPU N_VGetDeviceArrayPointer_Hip
#define N_VCopyFromDevice_GPU        N_VCopyFromDevice_Hip
#define GPU_DEVICE_SYNCHRONIZE       hipDeviceSynchronize
#define GPU_GET_LAST_ERROR           hipGetLastError
#define GPU_SUCCESS                  hipSuccess
#define GPU_GET_ERROR_NAME           hipGetErrorName
#endif

using GkoMatrixTypeGPU    = gko::matrix::Csr<sunrealtype, sunindextype>;
using SUNGkoMatrixTypeGPU = sundials::ginkgo::Matrix<GkoMatrixTypeGPU>;

N_Vector create_vector_gpu(sunindextype length, SUNContext sunctx)
{
  return N_VNew_GPU(length, sunctx);
}

__global__ void solution_kernel(const sunindextype nx, const sunindextype ny,
                                const sunrealtype dx, const sunrealtype dy,
                                const sunrealtype cos_sqr_t, sunrealtype* uarray)
{
  const sunindextype i = blockIdx.x * blockDim.x + threadIdx.x;
  const sunindextype j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
  {
    auto x = i * dx;
    auto y = j * dy;

    auto sin_sqr_x = sin(PI * x) * sin(PI * x);
    auto sin_sqr_y = sin(PI * y) * sin(PI * y);

    auto idx    = i + j * nx;
    uarray[idx] = sin_sqr_x * sin_sqr_y * cos_sqr_t + ONE;
  }
}

int Solution_gpu(sunrealtype t, N_Vector u, UserData& udata)
{
  const auto nx = udata.nx;
  const auto ny = udata.ny;
  const auto dx = udata.dx;
  const auto dy = udata.dy;

  const auto cos_sqr_t = cos(PI * t) * cos(PI * t);
  sunrealtype* uarray  = N_VGetDeviceArrayPointer_GPU(u);
  if (check_ptr(uarray, "N_VGetDeviceArrayPointer_GPU")) return -1;

  dim3 threads_per_block{16, 16};
  const auto nbx{(static_cast<unsigned int>(nx) + threads_per_block.x - 1) /
                 threads_per_block.x};
  const auto nby{(static_cast<unsigned int>(ny) + threads_per_block.y - 1) /
                 threads_per_block.y};
  dim3 num_blocks{nbx, nby};

  solution_kernel<<<num_blocks, threads_per_block>>>(nx, ny, dx, dy, cos_sqr_t,
                                                     uarray);

  GPU_DEVICE_SYNCHRONIZE();
  return GPU_GET_LAST_ERROR() == GPU_SUCCESS ? 0 : -1;
}

__global__ void f_kernel(const sunindextype nx, const sunindextype ny,
                         const sunrealtype dx, const sunrealtype dy,
                         const sunrealtype cx, const sunrealtype cy,
                         const sunrealtype cc, const sunrealtype bx,
                         const sunrealtype by, const sunrealtype sin_t_cos_t,
                         const sunrealtype cos_sqr_t, sunrealtype* uarray,
                         sunrealtype* farray)
{
  const sunindextype i = blockIdx.x * blockDim.x + threadIdx.x;
  const sunindextype j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
  {
    auto x = i * dx;
    auto y = j * dy;

    auto sin_sqr_x = sin(PI * x) * sin(PI * x);
    auto sin_sqr_y = sin(PI * y) * sin(PI * y);

    auto cos_sqr_x = cos(PI * x) * cos(PI * x);
    auto cos_sqr_y = cos(PI * y) * cos(PI * y);

    auto idx_c = i + j * nx;
    auto idx_n = i + (j + 1) * nx;
    auto idx_s = i + (j - 1) * nx;
    auto idx_e = (i + 1) + j * nx;
    auto idx_w = (i - 1) + j * nx;

    farray[idx_c] = cc * uarray[idx_c] + cx * (uarray[idx_w] + uarray[idx_e]) +
                    cy * (uarray[idx_s] + uarray[idx_n]) -
                    TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
                    bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
                    by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
  }
}

int f_gpu(sunrealtype t, N_Vector u, N_Vector f, UserData* udata)
{
  const auto nx = udata->nx;
  const auto ny = udata->ny;
  const auto dx = udata->dx;
  const auto dy = udata->dy;
  const auto kx = udata->kx;
  const auto ky = udata->ky;

  const sunrealtype cx   = kx / (dx * dx);
  const sunrealtype cy   = ky / (dy * dy);
  const sunrealtype cc   = -TWO * (cx + cy);
  const auto bx          = kx * TWO * PI * PI;
  const auto by          = ky * TWO * PI * PI;
  const auto sin_t_cos_t = sin(PI * t) * cos(PI * t);
  const auto cos_sqr_t   = cos(PI * t) * cos(PI * t);

  sunrealtype* uarray = N_VGetDeviceArrayPointer_GPU(u);
  if (check_ptr(uarray, "N_VGetDeviceArrayPointer_GPU")) return -1;
  sunrealtype* farray = N_VGetDeviceArrayPointer_GPU(f);
  if (check_ptr(farray, "N_VGetDeviceArrayPointer_GPU")) return -1;

  dim3 threads_per_block{16, 16};
  const auto nbx{(static_cast<unsigned int>(nx) + threads_per_block.x - 1) /
                 threads_per_block.x};
  const auto nby{(static_cast<unsigned int>(ny) + threads_per_block.y - 1) /
                 threads_per_block.y};
  dim3 num_blocks{nbx, nby};

  f_kernel<<<num_blocks, threads_per_block>>>(nx, ny, dx, dy, cx, cy, cc, bx,
                                              by, sin_t_cos_t, cos_sqr_t,
                                              uarray, farray);
  GPU_DEVICE_SYNCHRONIZE();
  return GPU_GET_LAST_ERROR() == GPU_SUCCESS ? 0 : -1;
}

__global__ void J_sn_kernel(const sunindextype nx, const sunindextype ny,
                            sunindextype* row_ptrs, sunindextype* col_idxs,
                            sunrealtype* mat_data)
{
  const sunindextype i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= 0 && i < nx)
  {
    mat_data[i] = ZERO;
    col_idxs[i] = i;
    row_ptrs[i] = i;

    auto col      = i + (ny - 1) * nx;
    auto idx      = (5 * (nx - 2) + 2) * (ny - 2) + nx + i;
    mat_data[idx] = ZERO;
    col_idxs[idx] = col;
    row_ptrs[col] = idx;
  }

  if (i == nx - 1) row_ptrs[nx * ny] = (5 * (nx - 2) + 2) * (ny - 2) + 2 * nx;
}

__global__ void J_we_kernel(const sunindextype nx, const sunindextype ny,
                            sunindextype* row_ptrs, sunindextype* col_idxs,
                            sunrealtype* mat_data)
{
  const sunindextype j = blockIdx.x * blockDim.x + threadIdx.x;

  if (j > 0 && j < ny - 1)
  {
    auto col      = j * nx;
    auto idx      = (5 * (nx - 2) + 2) * (j - 1) + nx;
    mat_data[idx] = ZERO;
    col_idxs[idx] = col;
    row_ptrs[col] = idx;

    col           = (nx - 1) + j * nx;
    idx           = (5 * (nx - 2) + 2) * (j - 1) + nx + 1 + 5 * (nx - 2);
    mat_data[idx] = ZERO;
    col_idxs[idx] = col;
    row_ptrs[col] = idx;
  }
}

__global__ void J_kernel(const sunindextype nx, const sunindextype ny,
                         const sunrealtype cx, const sunrealtype cy,
                         const sunrealtype cc, sunindextype* row_ptrs,
                         sunindextype* col_idxs, sunrealtype* mat_data)
{
  const sunindextype i = blockIdx.x * blockDim.x + threadIdx.x;
  const sunindextype j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
  {
    auto row   = i + j * nx;
    auto col_s = row - nx;
    auto col_w = row - 1;
    auto col_c = row;
    auto col_e = row + 1;
    auto col_n = row + nx;

    auto prior_nnz = (5 * (nx - 2) + 2) * (j - 1) + nx;
    auto idx       = prior_nnz + 1 + 5 * (i - 1);

    mat_data[idx]     = cy;
    mat_data[idx + 1] = cx;
    mat_data[idx + 2] = cc;
    mat_data[idx + 3] = cx;
    mat_data[idx + 4] = cy;

    col_idxs[idx]     = col_s;
    col_idxs[idx + 1] = col_w;
    col_idxs[idx + 2] = col_c;
    col_idxs[idx + 3] = col_e;
    col_idxs[idx + 4] = col_n;

    row_ptrs[row] = idx;
  }
}

int J_gpu(SUNMatrix J, UserData* udata)
{
  const auto nx = udata->nx;
  const auto ny = udata->ny;
  const auto dx = udata->dx;
  const auto dy = udata->dy;
  const auto kx = udata->kx;
  const auto ky = udata->ky;

  const sunrealtype cx = kx / (dx * dx);
  const sunrealtype cy = ky / (dy * dy);
  const sunrealtype cc = -TWO * (cx + cy);

  auto J_gko    = static_cast<SUNGkoMatrixTypeGPU*>(J->content)->GkoMtx();
  auto row_ptrs = J_gko->get_row_ptrs();
  auto col_idxs = J_gko->get_col_idxs();
  auto mat_data = J_gko->get_values();

  unsigned threads_per_block_bx = 16;
  unsigned num_blocks_bx = (nx + threads_per_block_bx - 1) / threads_per_block_bx;
  J_sn_kernel<<<num_blocks_bx, threads_per_block_bx>>>(nx, ny, row_ptrs,
                                                       col_idxs, mat_data);

  unsigned threads_per_block_by = 16;
  unsigned num_blocks_by = (ny + threads_per_block_by - 1) / threads_per_block_by;
  J_we_kernel<<<num_blocks_by, threads_per_block_by>>>(nx, ny, row_ptrs,
                                                       col_idxs, mat_data);

  dim3 threads_per_block_i{16, 16};
  const auto nbx{(static_cast<unsigned int>(nx) + threads_per_block_i.x - 1) /
                 threads_per_block_i.x};
  const auto nby{(static_cast<unsigned int>(ny) + threads_per_block_i.y - 1) /
                 threads_per_block_i.y};
  dim3 num_blocks_i{nbx, nby};
  J_kernel<<<num_blocks_i, threads_per_block_i>>>(nx, ny, cx, cy, cc, row_ptrs,
                                                  col_idxs, mat_data);

  GPU_DEVICE_SYNCHRONIZE();
  return GPU_GET_LAST_ERROR() == GPU_SUCCESS ? 0 : -1;
}

void copy_from_device_gpu(N_Vector v) { N_VCopyFromDevice_GPU(v); }

#undef N_VNew_GPU
#undef N_VGetDeviceArrayPointer_GPU
#undef N_VCopyFromDevice_GPU
#undef GPU_DEVICE_SYNCHRONIZE
#undef GPU_GET_LAST_ERROR
#undef GPU_SUCCESS
#undef GPU_GET_ERROR_NAME

// Include integrator, vector, matrix, and linear solver headers
#include <cvode/cvode.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
using GkoSolverType = gko::solver::Cg<sunrealtype>;

using SUNGkoMatrixType = sundials::ginkgo::Matrix<GkoMatrixType>;
using SUNGkoSolverType =
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>;

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data);

// Jacobian of RHS function
int J(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // SUNDIALS context
  sundials::Context sunctx;

  // ---------------
  // Setup UserData
  // ---------------

  // Read input options
  UserData udata;
  std::vector<std::string> args(argv + 1, argv + argc);
  if (ReadInputs(args, udata)) { return 1; }
  PrintUserData(udata);

  // ---------------------------------------
  // Create Ginkgo matrix and linear solver
  // ---------------------------------------

#if defined(USE_CUDA)
  auto gko_exec{gko::CudaExecutor::create(0, gko::OmpExecutor::create())};
#else
  auto gko_exec{gko::HipExecutor::create(0, gko::OmpExecutor::create())};
#endif

  udata.exec = gko_exec;

  // ---------------
  // Create vectors
  // ---------------

  // Create solution vector
  N_Vector u = create_vector_gpu(udata.nodes, sunctx);
  if (check_ptr(u, "N_VNew")) { return 1; }

  // Set initial condition
  int flag = Solution(ZERO, u, udata);
  if (check_flag(flag, "Solution")) { return 1; }

  // Create error vector
  N_Vector e = N_VClone(u);
  if (check_ptr(e, "N_VClone")) { return 1; }

  auto gko_matrix_dim = gko::dim<2>(udata.nodes, udata.nodes);
  auto gko_matrix_nnz{(5 * (udata.nx - 2) + 2) * (udata.ny - 2) + 2 * udata.nx};
  auto gko_matrix =
    gko::share(GkoMatrixType::create(gko_exec, gko_matrix_dim, gko_matrix_nnz));

  SUNGkoMatrixType A{gko_matrix, sunctx};

  // Use default stopping criteria
  auto crit{sundials::ginkgo::DefaultStop::build()
              .with_max_iters(static_cast<gko::uint64>(udata.liniters))
              .on(gko_exec)};

  // Use Jacobi preconditioner
  auto precon{
    gko::preconditioner::Jacobi<sunrealtype, sunindextype>::build().on(gko_exec)};

  auto gko_solver_factory = gko::share(GkoSolverType::build()
                                         .with_criteria(std::move(crit))
                                         .with_preconditioner(std::move(precon))
                                         .on(gko_exec));

  SUNGkoSolverType LS{gko_solver_factory, sunctx};

  // --------------
  // Setup CVODE
  // --------------

  // Create integrator
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }

  // Initialize integrator
  flag = CVodeInit(cvode_mem, f, ZERO, u);
  if (check_flag(flag, "CVodeInit")) { return 1; }

  // Specify tolerances
  flag = CVodeSStolerances(cvode_mem, udata.rtol, udata.atol);
  if (check_flag(flag, "CVodeSStolerances")) { return 1; }

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, (void*)&udata);
  if (check_flag(flag, "CVodeSetUserData")) { return 1; }

  // Attach linear solver
  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(flag, "CVodeSetLinearSolver")) { return 1; }

  // Attach Jacobian function
  flag = CVodeSetJacFn(cvode_mem, J);
  if (check_flag(flag, "CVodeSetJacFn")) { return 1; }

  // Set linear solver tolerance factor
  flag = CVodeSetEpsLin(cvode_mem, udata.epslin);
  if (check_flag(flag, "CVodeSetEpsLin")) { return 1; }

  // Set max steps between outputs
  flag = CVodeSetMaxNumSteps(cvode_mem, udata.maxsteps);
  if (check_flag(flag, "CVodeSetMaxNumSteps")) { return 1; }

  // Set stopping time
  flag = CVodeSetStopTime(cvode_mem, udata.tf);
  if (check_flag(flag, "CVodeSetStopTime")) { return 1; }

  // -----------------------
  // Loop over output times
  // -----------------------

  sunrealtype t     = ZERO;
  sunrealtype dTout = udata.tf / udata.nout;
  sunrealtype tout  = dTout;

  // Initial output
  flag = OpenOutput(udata);
  if (check_flag(flag, "OpenOutput")) { return 1; }

  flag = WriteOutput(t, u, e, udata);
  if (check_flag(flag, "WriteOutput")) { return 1; }

  for (int iout = 0; iout < udata.nout; iout++)
  {
    // Evolve in time
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_flag(flag, "CVode")) { break; }

    // Output solution and error
    flag = WriteOutput(t, u, e, udata);
    if (check_flag(flag, "WriteOutput")) { return 1; }

    // Update output time
    tout += dTout;
    tout = (tout > udata.tf) ? udata.tf : tout;
  }

  // Close output
  flag = CloseOutput(udata);
  if (check_flag(flag, "CloseOutput")) { return 1; }

  // --------------
  // Final outputs
  // --------------

  // Print final integrator stats
  std::cout << "Final integrator statistics:" << std::endl;
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) { return 1; }

  // Output final error
  flag = SolutionError(t, u, e, udata);
  if (check_flag(flag, "SolutionError")) { return 1; }

  sunrealtype maxerr = N_VMaxNorm(e);

  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<sunrealtype>::digits10)
            << "\nMax error = " << maxerr << std::endl;

  // --------------------
  // Clean up and return
  // --------------------

  udata.exec = nullptr;
  CVodeFree(&cvode_mem); // Free integrator memory
  N_VDestroy(u);         // Free vectors
  N_VDestroy(e);

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data)
{
  // Access problem data and set shortcuts
  auto udata = static_cast<UserData*>(user_data);

  // Initialize RHS vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  return f_gpu(t, u, f, udata);

  // Return success
  return 0;
}

// J routine to compute the ODE RHS Jacobian function df/dy(t,y). This
// explicitly set boundary entries to zero so J(t,y) has the same sparsity
// pattern as A = I - gamma * J(t,y).
int J(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  auto udata = static_cast<UserData*>(user_data);

  return J_gpu(J, udata);

  return 0;
}

//---- end of file ----
