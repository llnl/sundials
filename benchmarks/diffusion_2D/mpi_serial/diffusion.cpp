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
 * Serial diffusion functions
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

// Diffusion function
int laplacian(sunrealtype t, N_Vector u, N_Vector f, UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  int flag;
  sunindextype i, j;

  // Start exchange
  flag = udata->start_exchange(u);
  if (check_flag(&flag, "SendData", 1)) { return -1; }

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Constants for computing diffusion term
  sunrealtype cx = udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = -TWO * (cx + cy);

  // Access data arrays
  sunrealtype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  sunrealtype* farray = N_VGetArrayPointer(f);
  if (check_flag((void*)farray, "N_VGetArrayPointer", 0)) { return -1; }

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over subdomain and compute rhs forcing term
  if (udata->forcing)
  {
    sunrealtype x, y;
    sunrealtype sin_sqr_x, sin_sqr_y;
    sunrealtype cos_sqr_x, cos_sqr_y;

    sunrealtype bx = (udata->kx) * TWO * PI * PI;
    sunrealtype by = (udata->ky) * TWO * PI * PI;

    sunrealtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
    sunrealtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

    for (j = jstart; j < jend; j++)
    {
      for (i = istart; i < iend; i++)
      {
        x = (udata->is + i) * udata->dx;
        y = (udata->js + j) * udata->dy;

        sin_sqr_x = sin(PI * x) * sin(PI * x);
        sin_sqr_y = sin(PI * y) * sin(PI * y);

        cos_sqr_x = cos(PI * x) * cos(PI * x);
        cos_sqr_y = cos(PI * y) * cos(PI * y);

        farray[IDX(i, j, nx_loc)] =
          -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
          bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
          by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
      }
    }
  }

  // Iterate over subdomain interior and add rhs diffusion term
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + uarray[IDX(i, j + 1, nx_loc)]);
    }
  }

  // Wait for exchange receives
  flag = udata->end_exchange();
  if (check_flag(&flag, "UserData::end_excahnge", 1)) { return -1; }

  // Iterate over subdomain boundaries and add rhs diffusion term
  sunrealtype* Warray = udata->Wrecv;
  sunrealtype* Earray = udata->Erecv;
  sunrealtype* Sarray = udata->Srecv;
  sunrealtype* Narray = udata->Nrecv;

  // West face (updates south-west and north-west corners if necessary)
  if (udata->HaveNbrW)
  {
    i = 0;
    if (udata->HaveNbrS) // South-West corner
    {
      j = 0;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (Warray[j] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (Sarray[i] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (Warray[j] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    if (udata->HaveNbrN) // North-West corner
    {
      j = ny_loc - 1;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (Warray[j] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + Narray[i]);
    }
  }

  // East face (updates south-east and north-east corners if necessary)
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;
    if (udata->HaveNbrS) // South-East corner
    {
      j = 0;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + Earray[j]) +
        cy * (Sarray[i] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + Earray[j]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    if (udata->HaveNbrN) // North-East corner
    {
      j = ny_loc - 1;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + Earray[j]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + Narray[i]);
    }
  }

  // South face (excludes corners)
  if (udata->HaveNbrS)
  {
    j = 0;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (Sarray[i] + uarray[IDX(i, j + 1, nx_loc)]);
    }
  }

  // North face (excludes corners)
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + Narray[i]);
    }
  }

  // Return success
  return 0;
}

#if defined(USE_SUPERLU_DIST)
// Compute the global index of a node
//   i -- local x index
//   j -- local y index
//   x -- x processor coordinate
//   y -- y processor coordinate
static sunindextype global_index(sunindextype i, sunindextype j, int x, int y,
                                 UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Unpack values (same for all ranks)
  sunindextype qx = udata->qx;
  sunindextype qy = udata->qy;
  sunindextype rx = udata->rx;
  sunindextype ry = udata->ry;
  sunindextype ny = udata->ny;

  // offset from previous process rows
  sunindextype offset_p = ny *
                          ((qx + 1) * std::min((sunindextype)x, rx) +
                           qx * std::max((sunindextype)x - rx, (sunindextype)0));

  // offset within current process row
  sunindextype offset_c;
  if (x < rx)
  {
    offset_c = (qx + 1) * ((qy + 1) * std::min((sunindextype)y, ry) +
                           qy * std::max((sunindextype)y - ry, (sunindextype)0));
  }
  else
  {
    offset_c = qx * ((qy + 1) * std::min((sunindextype)y, ry) +
                     qy * std::max((sunindextype)y - ry, (sunindextype)0));
  }

  // nx_loc for the input (x,y) process (not necessarily the same as nx_loc in
  // user data because this could be the index for a neighboring rank)
  sunindextype nx_loc;
  if (x < rx) nx_loc = qx + 1;
  else nx_loc = qx;

  // global index for the requested local node
  return offset_p + offset_c + j * nx_loc + i;
}

// Compute a row in the global matrix
//   i -- local x index
//   j -- local y index
//   x -- x processor coordinate
//   y -- y processor coordinate
#if defined(BENCHMARK_ODE)
static int matrix_row(sunindextype i, sunindextype j, int x, int y,
                      UserData* udata, sunrealtype* vals, sunindextype* col_idx,
                      sunindextype* row_nnz)
#else
static int matrix_row(sunindextype i, sunindextype j, int x, int y,
                      UserData* udata, sunrealtype cj, sunrealtype* vals,
                      sunindextype* col_idx, sunindextype* row_nnz)
#endif
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Unpack values
  int npx = udata->npx;
  int npy = udata->npy;

  sunindextype qx = udata->qx;
  sunindextype qy = udata->qy;
  sunindextype rx = udata->rx;
  sunindextype ry = udata->ry;

  // -----------------
  // global boundaries
  // -----------------

  if ((y == 0 && j == 0) ||                       // south
      (x == 0 && i == 0) ||                       // west
      (x == npx - 1 && i == udata->nx_loc - 1) || // east
      (y == npy - 1 && j == udata->ny_loc - 1))   // north
  {
#if defined(BENCHMARK_ODE)
    vals[0] = 0;
#else
    vals[0] = cj;
#endif
    col_idx[0] = global_index(i, j, x, y, udata);
    *row_nnz   = 1;
    return 0;
  }

  // --------
  // interior
  // --------

  // Ordering of columns for output
  //   0. West Neighbor
  //   1. South Neighbor
  //   Interior:
  //      2. South
  //      3. West
  //      4. Center
  //      5. East
  //      6. North
  //   7. North Neighbor
  //   8. East Neighbor

  // Constants for computing Jacobian
#if defined(BENCHMARK_ODE)
  sunrealtype cx = udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = -TWO * (cx + cy);
#else
  sunrealtype cx = -udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = -udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = cj - TWO * (cx + cy);
#endif

  // Value and Columns index
  sunindextype idx = 0;

  // west neighbor
  if (i == 0)
  {
    // neighbor nx_loc
    sunindextype nx_loc;
    nx_loc = (x - 1 < rx) ? qx + 1 : qx;

    vals[idx]    = cx;
    col_idx[idx] = global_index(nx_loc - 1, j, x - 1, y, udata);
    idx += 1;
  }

  // south neighbor
  if (j == 0)
  {
    // neighbor ny_loc
    sunindextype ny_loc;
    ny_loc = (y - 1 < ry) ? qy + 1 : qy;

    vals[idx]    = cy;
    col_idx[idx] = global_index(i, ny_loc - 1, x, y - 1, udata);
    idx += 1;
  }

  // south interior
  if (j > 0)
  {
    vals[idx]    = cy;
    col_idx[idx] = global_index(i, j - 1, x, y, udata);
    idx += 1;
  }

  // west interior
  if (i > 0)
  {
    vals[idx]    = cx;
    col_idx[idx] = global_index(i - 1, j, x, y, udata);
    idx += 1;
  }

  // center
  vals[idx]    = cc;
  col_idx[idx] = global_index(i, j, x, y, udata);
  idx += 1;

  // east interior
  if (i < udata->nx_loc - 1)
  {
    vals[idx]    = cx;
    col_idx[idx] = global_index(i + 1, j, x, y, udata);
    idx += 1;
  }

  // north interior
  if (j < udata->ny_loc - 1)
  {
    vals[idx]    = cy;
    col_idx[idx] = global_index(i, j + 1, x, y, udata);
    idx += 1;
  }

  // north neighbor
  if (j == udata->ny_loc - 1)
  {
    vals[idx]    = cy;
    col_idx[idx] = global_index(i, 0, x, y + 1, udata);
    idx += 1;
  }

  // east neighbor
  if (i == udata->nx_loc - 1)
  {
    vals[idx]    = cx;
    col_idx[idx] = global_index(0, j, x + 1, y, udata);
    idx += 1;
  }

  // non-zeros in this row
  *row_nnz = 5;

  return 0;
}

#if defined(BENCHMARK_ODE)
int laplacian_matrix_sludist(N_Vector u, SUNMatrix L, UserData* udata)
#else
int laplacian_matrix_sludist(N_Vector u, sunrealtype cj, SUNMatrix L,
                             UserData* udata)
#endif
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Set shortcuts
  SuperMatrix* Lsuper    = SUNMatrix_SLUNRloc_SuperMatrix(L);
  NRformat_loc* Lstore   = (NRformat_loc*)Lsuper->Store;
  sunindextype* row_ptrs = Lstore->rowptr;
  sunindextype* col_idxs = Lstore->colind;
  sunrealtype* data      = (sunrealtype*)Lstore->nzval;

  int x = udata->idx;
  int y = udata->idy;

  // Initial row pointer value
  row_ptrs[0] = 0;

  sunindextype idx     = 0;
  sunindextype row_nnz = 0;

  for (sunindextype j = 0; j < udata->ny_loc; j++)
  {
    for (sunindextype i = 0; i < udata->nx_loc; i++)
    {
#if defined(BENCHMARK_ODE)
      matrix_row(i, j, x, y, udata, data + row_ptrs[idx],
                 col_idxs + row_ptrs[idx], &row_nnz);
#else
      matrix_row(i, j, x, y, udata, cj, data + row_ptrs[idx],
                 col_idxs + row_ptrs[idx], &row_nnz);
#endif
      row_ptrs[idx + 1] = row_ptrs[idx] + row_nnz;
      idx++;
    }
  }

  Lstore->nnz_loc = row_ptrs[idx];
  Lstore->m_loc   = udata->nx_loc * udata->ny_loc;
  Lstore->fst_row = global_index(0, 0, x, y, udata);

  // Return success
  return 0;
}
#endif
