/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                Cody J. Balos @ LLNL
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
 * A simple implementation of a parallel structured Cartesian mesh class that
 * supports up to 3 spatial dimensions and an arbitrary number of degrees of
 * freedom, and that uses Kokkos views to store communication buffer data.
 * ----------------------------------------------------------------------------*/

#ifndef _KOKKOSPARGRID_H
#define _KOKKOSPARGRID_H

#include <Kokkos_Core.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sundials/sundials_types.h>

/* Set Kokkos execution space and type shortcuts */
#if defined(USE_CUDA)
using ExecSpace = Kokkos::Cuda;
using MemSpace  = Kokkos::CudaSpace;
#elif defined(USE_HIP)
#if KOKKOS_VERSION / 10000 > 3
using ExecSpace = Kokkos::HIP;
using MemSpace  = Kokkos::HIPSpace;
#else
using ExecSpace = Kokkos::Experimental::HIP;
using MemSpace  = Kokkos::Experimental::HIPSpace;
#endif
#elif defined(USE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
using MemSpace  = Kokkos::HostSpace;
#else
using ExecSpace = Kokkos::Serial;
using MemSpace  = Kokkos::HostSpace;
#endif
using Vec1D     = Kokkos::View<sunrealtype*, MemSpace>;
using Vec4D     = Kokkos::View<sunrealtype****, MemSpace>;
using Vec1DHost = Vec1D::HostMirror;
using Vec4DHost = Vec4D::HostMirror;
using Range3D   = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>;

namespace sundials_tools {

// Types of boundaries supported.
enum class BoundaryType
{
  PERIODIC
};

// Types of stencils supported.
enum class StencilType
{
  UPWIND
};

template<typename GLOBALINT>
class ParallelGrid
{
public:
  // Constructor that creates a new ParallelGrid object.
  // [in] - the memory helper to use for allocating the MPI buffers
  // [in,out] comm - on input, the overall MPI communicator, on output, the cartesian communicator
  // [in] a[] - an array of length 3 which defines the domain [a,b]
  // [in] b[] - an array of length 3 which defines the domain [a,b]
  // [in] npts[] - an array of length 3 which defines the number of mesh points in each dimension
  // [in] dof - the number of degrees of freedom in each dimension
  // [in] bc - the type of boundary conditions (see BoundaryType)
  // [in] st - the stencil to use (see StencilType)
  // [in] npxyz - the number of processors in each dimension; defaults to 0 which means MPI will choose
  // [in] reorder - should MPI_Cart_create do process reordering to optimize or not; defaults to false (some MPI implementations ignore this)
  ParallelGrid(MPI_Comm* comm, const sunrealtype a[], const sunrealtype b[],
               const GLOBALINT npts[], int dof_, BoundaryType bc_,
               StencilType st_, const sunrealtype c,
               const int npxyz[] = nullptr, bool reorder = false)
    : nx(1),
      ny(1),
      nz(1),
      nxl(1),
      nyl(1),
      nzl(1),
      npx(1),
      npy(1),
      npz(1),
      dx(0.0),
      dy(0.0),
      dz(0.0),
      ax(0.0),
      ay(0.0),
      az(0.0),
      bx(0.0),
      by(0.0),
      bz(0.0),
      dof(dof_),
      upwindRight(true),
      dims{0, 0, 0},
      coords{0, 0, 0},
      bc(bc_),
      st(st_)
  {
    assert(st == StencilType::UPWIND);

    /* Set up MPI Cartesian communicator */
    if (npxyz)
    {
      dims[0] = npxyz[0];
      dims[1] = npxyz[1];
      dims[2] = npxyz[2];
    }

    int retval, nprocs;
#ifdef NDEBUG
    // Suppress unused variable warning
    ((void)retval);
#endif
    MPI_Comm_size(*comm, &nprocs);
    retval = MPI_Dims_create(nprocs, 3, dims);
    assert(retval == MPI_SUCCESS);

    int periods[] = {bc == BoundaryType::PERIODIC, bc == BoundaryType::PERIODIC,
                     bc == BoundaryType::PERIODIC};
    retval        = MPI_Cart_create(*comm, 3, dims, periods, reorder, comm);
    assert(retval == MPI_SUCCESS);

    retval = MPI_Cart_get(*comm, 3, dims, periods, coords);
    assert(retval == MPI_SUCCESS);

    cart_comm = *comm;

    /* Set upwinding direction */
    upwindRight = (c > 0.0);

    /* Set up information for the first spatial dimension */
    npx    = dims[0];
    nx     = npts[0];
    ax     = a[0];
    bx     = b[0];
    dx     = (bx - ax) / (sunrealtype)nx;
    int is = nx * (coords[0]) / npx;
    int ie = nx * (coords[0] + 1) / npx - 1;
    nxl    = ie - is + 1;
    neq    = dof * nxl;

    /* Set up information for the second spatial dimension */
    npy    = dims[1];
    ny     = npts[1];
    ay     = a[1];
    by     = b[1];
    dy     = (by - ay) / (sunrealtype)ny;
    int js = ny * (coords[1]) / npy;
    int je = ny * (coords[1] + 1) / npy - 1;
    nyl    = je - js + 1;
    neq *= nyl;

    /* Set up information for the third spatial dimension */
    npz    = dims[2];
    nz     = npts[2];
    az     = a[2];
    bz     = b[2];
    dz     = (bz - az) / (sunrealtype)nz;
    int ks = nz * (coords[2]) / npz;
    int ke = nz * (coords[2] + 1) / npz - 1;
    nzl    = ke - ks + 1;
    neq *= nzl;

    /* Allocate buffers for nearest-neighbor exchange */
    if (st == StencilType::UPWIND) { AllocateBuffersUpwind(); }
  }

  // TODO:
  //  - support non-periodic boundary conditions
  // For all faces where neighbors exist: determine neighbor process indices.
  // For all faces: allocate upwind exchange buffers.
  void AllocateBuffersUpwind()
  {
    int retval;
#ifdef NDEBUG
    // Suppress unused variable warning
    ((void)retval);
#endif

    /* Allocate send/receive buffers and determine ID for communication West */
    if (upwindRight)
    {
      Wrecv_  = Vec1D("Wrecv", dof * nyl * nzl);
      WrecvH_ = Kokkos::create_mirror_view(Wrecv_);
    }
    else
    {
      Wsend_  = Vec1D("Wsend", dof * nyl * nzl);
      WsendH_ = Kokkos::create_mirror_view(Wsend_);
    }
    ipW = MPI_PROC_NULL;
    if ((coords[0] > 0) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0] - 1, coords[1], coords[2]};
      retval         = MPI_Cart_rank(cart_comm, nbcoords, &ipW);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication East */
    if (upwindRight)
    {
      Esend_  = Vec1D("Esend", dof * nyl * nzl);
      EsendH_ = Kokkos::create_mirror_view(Esend_);
    }
    else
    {
      Erecv_  = Vec1D("Erecv", dof * nyl * nzl);
      ErecvH_ = Kokkos::create_mirror_view(Erecv_);
    }
    ipE = MPI_PROC_NULL;
    if ((coords[0] < dims[0] - 1) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0] + 1, coords[1], coords[2]};
      retval         = MPI_Cart_rank(cart_comm, nbcoords, &ipE);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication South */
    if (upwindRight)
    {
      Srecv_  = Vec1D("Srecv", dof * nxl * nzl);
      SrecvH_ = Kokkos::create_mirror_view(Srecv_);
    }
    else
    {
      Ssend_  = Vec1D("Ssend", dof * nxl * nzl);
      SsendH_ = Kokkos::create_mirror_view(Ssend_);
    }
    ipS = MPI_PROC_NULL;
    if ((coords[1] > 0) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1] - 1, coords[2]};
      retval         = MPI_Cart_rank(cart_comm, nbcoords, &ipS);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication North */
    if (upwindRight)
    {
      Nsend_  = Vec1D("Nsend", dof * nxl * nzl);
      NsendH_ = Kokkos::create_mirror_view(Nsend_);
    }
    else
    {
      Nrecv_  = Vec1D("Nrecv", dof * nxl * nzl);
      NrecvH_ = Kokkos::create_mirror_view(Nrecv_);
    }
    ipN = MPI_PROC_NULL;
    if ((coords[1] < dims[1] - 1) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1] + 1, coords[2]};
      retval         = MPI_Cart_rank(cart_comm, nbcoords, &ipN);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication Back */
    if (upwindRight)
    {
      Brecv_  = Vec1D("Brecv", dof * nxl * nyl);
      BrecvH_ = Kokkos::create_mirror_view(Brecv_);
    }
    else
    {
      Bsend_  = Vec1D("Bsend", dof * nxl * nyl);
      BsendH_ = Kokkos::create_mirror_view(Bsend_);
    }
    ipB = MPI_PROC_NULL;
    if ((coords[2] > 0) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1], coords[2] - 1};
      retval         = MPI_Cart_rank(cart_comm, nbcoords, &ipB);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication Front */
    if (upwindRight)
    {
      Fsend_  = Vec1D("Fsend", dof * nxl * nyl);
      FsendH_ = Kokkos::create_mirror_view(Fsend_);
    }
    else
    {
      Frecv_  = Vec1D("Frecv", dof * nxl * nyl);
      FrecvH_ = Kokkos::create_mirror_view(Frecv_);
    }
    ipF = MPI_PROC_NULL;
    if ((coords[2] < dims[2] - 1) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1], coords[2] + 1};
      retval         = MPI_Cart_rank(cart_comm, nbcoords, &ipF);
      assert(retval == MPI_SUCCESS);
    }
  }

  // Initiate non-blocking neighbor communication
  int ExchangeStart()
  {
    int retval = 0;
    nreq       = 0;

    // Initialize all requests in array
    for (int i = 0; i < 12; i++) { req[i] = MPI_REQUEST_NULL; }

    // Open an Irecv buffer on host for each neighbor
    if ((ipW != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(WrecvH_.data(), dof * nyl * nzl, MPI_SUNREALTYPE, ipW,
                         1, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipE != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(ErecvH_.data(), dof * nyl * nzl, MPI_SUNREALTYPE, ipE,
                         0, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipS != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(SrecvH_.data(), dof * nxl * nzl, MPI_SUNREALTYPE, ipS,
                         3, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipN != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(NrecvH_.data(), dof * nxl * nzl, MPI_SUNREALTYPE, ipN,
                         2, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipB != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(BrecvH_.data(), dof * nxl * nyl, MPI_SUNREALTYPE, ipB,
                         5, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipF != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(FrecvH_.data(), dof * nxl * nyl, MPI_SUNREALTYPE, ipF,
                         4, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    // Send data to neighbors, first copying from device to host buffers
    if ((ipW != MPI_PROC_NULL) && (!upwindRight))
    {
      Kokkos::deep_copy(WsendH_, Wsend_);
      retval = MPI_Isend(WsendH_.data(), dof * nyl * nzl, MPI_SUNREALTYPE, ipW,
                         0, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipE != MPI_PROC_NULL) && (upwindRight))
    {
      Kokkos::deep_copy(EsendH_, Esend_);
      retval = MPI_Isend(EsendH_.data(), dof * nyl * nzl, MPI_SUNREALTYPE, ipE,
                         1, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipS != MPI_PROC_NULL) && (!upwindRight))
    {
      Kokkos::deep_copy(SsendH_, Ssend_);
      retval = MPI_Isend(SsendH_.data(), dof * nxl * nzl, MPI_SUNREALTYPE, ipS,
                         2, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipN != MPI_PROC_NULL) && (upwindRight))
    {
      Kokkos::deep_copy(NsendH_, Nsend_);
      retval = MPI_Isend(NsendH_.data(), dof * nxl * nzl, MPI_SUNREALTYPE, ipN,
                         3, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipB != MPI_PROC_NULL) && (!upwindRight))
    {
      Kokkos::deep_copy(BsendH_, Bsend_);
      retval = MPI_Isend(BsendH_.data(), dof * nxl * nyl, MPI_SUNREALTYPE, ipB,
                         4, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipF != MPI_PROC_NULL) && (upwindRight))
    {
      Kokkos::deep_copy(FsendH_, Fsend_);
      retval = MPI_Isend(FsendH_.data(), dof * nxl * nyl, MPI_SUNREALTYPE, ipF,
                         5, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    return retval;
  }

  // Waits for neighbor exchange to finish.
  int ExchangeEnd()
  {
    MPI_Status stat[12];
    int retval;

    // return automatically with success if there are no outstanding requests
    if (nreq == 0) { return (0); }

    // Wait for messages to finish send/receive
    retval = MPI_Waitall(nreq, req, stat);
    assert(retval == MPI_SUCCESS);

    // Copy data from host to device buffers
    if ((ipW != MPI_PROC_NULL) && (upwindRight))
      Kokkos::deep_copy(Wrecv_, WrecvH_);
    if ((ipE != MPI_PROC_NULL) && (!upwindRight))
      Kokkos::deep_copy(Erecv_, ErecvH_);
    if ((ipS != MPI_PROC_NULL) && (upwindRight))
      Kokkos::deep_copy(Srecv_, SrecvH_);
    if ((ipN != MPI_PROC_NULL) && (!upwindRight))
      Kokkos::deep_copy(Nrecv_, NrecvH_);
    if ((ipB != MPI_PROC_NULL) && (upwindRight))
      Kokkos::deep_copy(Brecv_, BrecvH_);
    if ((ipF != MPI_PROC_NULL) && (!upwindRight))
      Kokkos::deep_copy(Frecv_, FrecvH_);

    return retval;
  }

  // Prints out information about the ParallelGrid to stdout.
  void PrintInfo()
  {
    printf("ParallelGrid Info:\n");
    printf("    dimensions = %d\n", 3);
    printf("    processors = {%d, %d, %d}\n", npx, npy, npz);
    printf("        domain = {[%g,%g], [%g,%g], [%g,%g]}\n", ax, bx, ay, by, az,
           bz);
    printf("   global npts = {%li, %li, %li}\n", (long int)nx, (long int)ny,
           (long int)nz);
    printf("    local npts = {%d, %d, %d}\n", nxl, nyl, nzl);
    printf("  mesh spacing = {%g, %g, %g}\n", dx, dy, dz);
    if (upwindRight) { printf("    upwind dir = right\n"); }
    else { printf("    upwind dir = left\n"); }
  }

  // Saves the mesh to a file.
  //    First row is x. Second row is y. Third row is z.
  //    Can be loaded into MATLAB like so:
  //      mesh = loadtxt('mesh.txt');
  //      [X,Y,Z] = meshgrid(mesh(1,:),mesh(2,:),mesh(3,:));
  void MeshToFile(const std::string& fname)
  {
    std::ofstream mesh_file;
    mesh_file.open(fname);
    mesh_file << std::setprecision(16);
    for (GLOBALINT i = 0; i < nx; i++) { mesh_file << " " << dx * i; }
    mesh_file << std::endl;
    for (GLOBALINT i = 0; i < ny; i++) { mesh_file << " " << dy * i; }
    mesh_file << std::endl;
    for (GLOBALINT i = 0; i < nz; i++) { mesh_file << " " << dz * i; }
    mesh_file << std::endl;
    mesh_file.close();
  }

  int nprocs() const { return npx * npy * npz; }

  GLOBALINT npts() const { return nx * ny * nz; }

  GLOBALINT nptsl() const { return nxl * nyl * nzl; }

  GLOBALINT neql() const { return dof * nptsl(); }

  sunrealtype* GetRecvView(const std::string& direction)
  {
    if (direction == "WEST")
    {
      return static_cast<sunrealtype*>(Wrecv_.data());
    }
    else if (direction == "EAST")
    {
      return static_cast<sunrealtype*>(Erecv_.data());
    }
    else if (direction == "NORTH")
    {
      return static_cast<sunrealtype*>(Nrecv_.data());
    }
    else if (direction == "SOUTH")
    {
      return static_cast<sunrealtype*>(Srecv_.data());
    }
    else if (direction == "FRONT")
    {
      return static_cast<sunrealtype*>(Frecv_.data());
    }
    else if (direction == "BACK")
    {
      return static_cast<sunrealtype*>(Brecv_.data());
    }
    else
    {
      assert(direction == "ILLEGAL");
      return nullptr;
    }
  }

  sunrealtype* GetSendView(const std::string& direction)
  {
    if (direction == "WEST")
    {
      return static_cast<sunrealtype*>(Wsend_.data());
    }
    else if (direction == "EAST")
    {
      return static_cast<sunrealtype*>(Esend_.data());
    }
    else if (direction == "NORTH")
    {
      return static_cast<sunrealtype*>(Nsend_.data());
    }
    else if (direction == "SOUTH")
    {
      return static_cast<sunrealtype*>(Ssend_.data());
    }
    else if (direction == "FRONT")
    {
      return static_cast<sunrealtype*>(Fsend_.data());
    }
    else if (direction == "BACK")
    {
      return static_cast<sunrealtype*>(Bsend_.data());
    }
    else
    {
      assert(direction == "ILLEGAL");
      return nullptr;
    }
  }

  GLOBALINT nx, ny, nz;   /* number of intervals globally       */
  int nxl, nyl, nzl;      /* number of intervals locally        */
  int npx, npy, npz;      /* number of processes                */
  sunrealtype dx, dy, dz; /* mesh spacing                       */
  sunrealtype ax, ay, az; /* domain in [a, b]                   */
  sunrealtype bx, by, bz;
  int dof; /* degrees of freedom per node        */
  int neq; /* total number of equations locally  */

  int ipW, ipE; /* MPI ranks for neighbor procs       */
  int ipS, ipN;
  int ipB, ipF;
  bool upwindRight; /* Upwind dir: true/false == R/L      */

  int dims[3];
  int coords[3];

private:
  MPI_Comm cart_comm; /* MPI cartesian communicator         */
  MPI_Request req[12];
  int nreq;

  BoundaryType bc;
  StencilType st;

  Vec1D Wsend_; /* MPI send/recv buffers              */
  Vec1D Esend_;
  Vec1D Ssend_;
  Vec1D Nsend_;
  Vec1D Bsend_;
  Vec1D Fsend_;
  Vec1D Wrecv_;
  Vec1D Erecv_;
  Vec1D Srecv_;
  Vec1D Nrecv_;
  Vec1D Brecv_;
  Vec1D Frecv_;
  Vec1DHost WsendH_; /* MPI send/recv buffers (host)       */
  Vec1DHost EsendH_;
  Vec1DHost SsendH_;
  Vec1DHost NsendH_;
  Vec1DHost BsendH_;
  Vec1DHost FsendH_;
  Vec1DHost WrecvH_;
  Vec1DHost ErecvH_;
  Vec1DHost SrecvH_;
  Vec1DHost NrecvH_;
  Vec1DHost BrecvH_;
  Vec1DHost FrecvH_;
};

} // namespace sundials_tools

#endif
