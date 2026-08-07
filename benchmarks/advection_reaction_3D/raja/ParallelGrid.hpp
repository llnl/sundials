/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * A simple implementation of a parallel structured cartesian mesh class that
 * supports up to 3 dimensions and an arbitrary number of degrees of freedom.
 * ----------------------------------------------------------------------------*/

#ifndef _SIMPLEPARGRID_H
#define _SIMPLEPARGRID_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sundials/sundials_memory.h>

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

template<typename REAL, typename GLOBALINT>
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
  // [in] width - the stencil width; defaults to 1
  // [in] npxyz - the number of processors in each dimension; defaults to 0 which means MPI will choose
  // [in] reorder - should MPI_Cart_create do process reordering to optimize or not; defaults to false (some MPI implementations ignore this)
  ParallelGrid(SUNMemoryHelper memhelp, MPI_Comm* comm, const REAL a[],
               const REAL b[], const GLOBALINT npts[], int dof, BoundaryType bc,
               StencilType st, const REAL c, int width = 1,
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
      dof(dof),
      dims{0, 0, 0},
      coords{0, 0, 0},
      bc(bc),
      st(st),
      width(width),
      upwindRight(true),
      memhelp(memhelp)

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
    dx     = (bx - ax) / (REAL)nx;
    int is = nx * (coords[0]) / npx;
    int ie = nx * (coords[0] + 1) / npx - 1;
    nxl    = ie - is + 1;
    neq    = dof * nxl;

    /* Set up information for the second spatial dimension */
    npy    = dims[1];
    ny     = npts[1];
    ay     = a[1];
    by     = b[1];
    dy     = (by - ay) / (REAL)ny;
    int js = ny * (coords[1]) / npy;
    int je = ny * (coords[1] + 1) / npy - 1;
    nyl    = je - js + 1;
    neq *= nyl;

    /* Set up information for the third spatial dimension */
    npz    = dims[2];
    nz     = npts[2];
    az     = a[2];
    bz     = b[2];
    dz     = (bz - az) / (REAL)nz;
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
  // For all faces: allocate exchange buffers.
  void AllocateBuffersUpwind()
  {
    /* Allocate send/receive buffers and determine ID for communication West */
    if (upwindRight)
    {
      SUNMemoryHelper_Alloc(memhelp, &Wrecv_,
                            sizeof(REAL) * dof * width * nyl * nzl,
                            memoryType(), nullptr);
    }
    else
    {
      SUNMemoryHelper_Alloc(memhelp, &Wsend_,
                            sizeof(REAL) * dof * width * nyl * nzl,
                            memoryType(), nullptr);
    }
    ipW = MPI_PROC_NULL;
    if ((coords[0] > 0) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0] - 1, coords[1], coords[2]};
      int retval     = MPI_Cart_rank(cart_comm, nbcoords, &ipW);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication East */
    if (upwindRight)
    {
      SUNMemoryHelper_Alloc(memhelp, &Esend_,
                            sizeof(REAL) * dof * width * nyl * nzl,
                            memoryType(), nullptr);
    }
    else
    {
      SUNMemoryHelper_Alloc(memhelp, &Erecv_,
                            sizeof(REAL) * dof * width * nyl * nzl,
                            memoryType(), nullptr);
    }
    ipE = MPI_PROC_NULL;
    if ((coords[0] < dims[0] - 1) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0] + 1, coords[1], coords[2]};
      int retval     = MPI_Cart_rank(cart_comm, nbcoords, &ipE);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication South */
    if (upwindRight)
    {
      SUNMemoryHelper_Alloc(memhelp, &Srecv_,
                            sizeof(REAL) * dof * width * nxl * nzl,
                            memoryType(), nullptr);
    }
    else
    {
      SUNMemoryHelper_Alloc(memhelp, &Ssend_,
                            sizeof(REAL) * dof * width * nxl * nzl,
                            memoryType(), nullptr);
    }
    ipS = MPI_PROC_NULL;
    if ((coords[1] > 0) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1] - 1, coords[2]};
      int retval     = MPI_Cart_rank(cart_comm, nbcoords, &ipS);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication North */
    if (upwindRight)
    {
      SUNMemoryHelper_Alloc(memhelp, &Nsend_,
                            sizeof(REAL) * dof * width * nxl * nzl,
                            memoryType(), nullptr);
    }
    else
    {
      SUNMemoryHelper_Alloc(memhelp, &Nrecv_,
                            sizeof(REAL) * dof * width * nxl * nzl,
                            memoryType(), nullptr);
    }
    ipN = MPI_PROC_NULL;
    if ((coords[1] < dims[1] - 1) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1] + 1, coords[2]};
      int retval     = MPI_Cart_rank(cart_comm, nbcoords, &ipN);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication Back */
    if (upwindRight)
    {
      SUNMemoryHelper_Alloc(memhelp, &Brecv_,
                            sizeof(REAL) * dof * width * nxl * nyl,
                            memoryType(), nullptr);
    }
    else
    {
      SUNMemoryHelper_Alloc(memhelp, &Bsend_,
                            sizeof(REAL) * dof * width * nxl * nyl,
                            memoryType(), nullptr);
    }
    ipB = MPI_PROC_NULL;
    if ((coords[2] > 0) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1], coords[2] - 1};
      int retval     = MPI_Cart_rank(cart_comm, nbcoords, &ipB);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication Front */
    if (upwindRight)
    {
      SUNMemoryHelper_Alloc(memhelp, &Fsend_,
                            sizeof(REAL) * dof * width * nxl * nyl,
                            memoryType(), nullptr);
    }
    else
    {
      SUNMemoryHelper_Alloc(memhelp, &Frecv_,
                            sizeof(REAL) * dof * width * nxl * nyl,
                            memoryType(), nullptr);
    }
    ipF = MPI_PROC_NULL;
    if ((coords[2] < dims[2] - 1) || (bc == BoundaryType::PERIODIC))
    {
      int nbcoords[] = {coords[0], coords[1], coords[2] + 1};
      int retval     = MPI_Cart_rank(cart_comm, nbcoords, &ipF);
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

    // Open an Irecv buffer for each neighbor
    if ((ipW != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(getRecvBuffer("WEST"), dof * nyl * nzl,
                         MPI_SUNREALTYPE, ipW, 1, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipE != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(getRecvBuffer("EAST"), dof * nyl * nzl,
                         MPI_SUNREALTYPE, ipE, 0, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipS != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(getRecvBuffer("SOUTH"), dof * nxl * nzl,
                         MPI_SUNREALTYPE, ipS, 3, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipN != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(getRecvBuffer("NORTH"), dof * nxl * nzl,
                         MPI_SUNREALTYPE, ipN, 2, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipB != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(getRecvBuffer("BACK"), dof * nxl * nyl,
                         MPI_SUNREALTYPE, ipB, 5, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipF != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(getRecvBuffer("FRONT"), dof * nxl * nyl,
                         MPI_SUNREALTYPE, ipF, 4, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    // Send data to neighbors
    if ((ipW != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Isend(getSendBuffer("WEST"), dof * nyl * nzl,
                         MPI_SUNREALTYPE, ipW, 0, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipE != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Isend(getSendBuffer("EAST"), dof * nyl * nzl,
                         MPI_SUNREALTYPE, ipE, 1, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipS != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Isend(getSendBuffer("SOUTH"), dof * nxl * nzl,
                         MPI_SUNREALTYPE, ipS, 2, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipN != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Isend(getSendBuffer("NORTH"), dof * nxl * nzl,
                         MPI_SUNREALTYPE, ipN, 3, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipB != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Isend(getSendBuffer("BACK"), dof * nxl * nyl,
                         MPI_SUNREALTYPE, ipB, 4, cart_comm, req + nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipF != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Isend(getSendBuffer("FRONT"), dof * nxl * nyl,
                         MPI_SUNREALTYPE, ipF, 5, cart_comm, req + nreq);
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

  REAL* getRecvBuffer(const std::string& direction)
  {
    if (direction == "WEST") { return static_cast<REAL*>(Wrecv_->ptr); }
    else if (direction == "EAST") { return static_cast<REAL*>(Erecv_->ptr); }
    else if (direction == "NORTH") { return static_cast<REAL*>(Nrecv_->ptr); }
    else if (direction == "SOUTH") { return static_cast<REAL*>(Srecv_->ptr); }
    else if (direction == "FRONT") { return static_cast<REAL*>(Frecv_->ptr); }
    else if (direction == "BACK") { return static_cast<REAL*>(Brecv_->ptr); }
    else
    {
      assert(direction == "ILLEGAL");
      return nullptr;
    }
  }

  REAL* getSendBuffer(const std::string& direction)
  {
    if (direction == "WEST") { return static_cast<REAL*>(Wsend_->ptr); }
    else if (direction == "EAST") { return static_cast<REAL*>(Esend_->ptr); }
    else if (direction == "NORTH") { return static_cast<REAL*>(Nsend_->ptr); }
    else if (direction == "SOUTH") { return static_cast<REAL*>(Ssend_->ptr); }
    else if (direction == "FRONT") { return static_cast<REAL*>(Fsend_->ptr); }
    else if (direction == "BACK") { return static_cast<REAL*>(Bsend_->ptr); }
    else
    {
      assert(direction == "ILLEGAL");
      return nullptr;
    }
  }

  ~ParallelGrid()
  {
    if (upwindRight)
    {
      SUNMemoryHelper_Dealloc(memhelp, Esend_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Nsend_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Fsend_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Wrecv_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Srecv_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Brecv_, nullptr);
    }
    else
    {
      SUNMemoryHelper_Dealloc(memhelp, Wsend_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Ssend_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Bsend_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Erecv_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Nrecv_, nullptr);
      SUNMemoryHelper_Dealloc(memhelp, Frecv_, nullptr);
    }
  }

  GLOBALINT nx, ny, nz; /* number of intervals globally       */
  int nxl, nyl, nzl;    /* number of intervals locally        */
  int npx, npy, npz;    /* number of processes                */
  REAL dx, dy, dz;      /* mesh spacing                       */
  REAL ax, ay, az;      /* domain in [a, b]                   */
  REAL bx, by, bz;
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
  int width;

  SUNMemoryHelper memhelp;
  SUNMemory Wsend_; /* MPI send/recv buffers              */
  SUNMemory Esend_;
  SUNMemory Ssend_;
  SUNMemory Nsend_;
  SUNMemory Bsend_;
  SUNMemory Fsend_;
  SUNMemory Wrecv_;
  SUNMemory Erecv_;
  SUNMemory Srecv_;
  SUNMemory Nrecv_;
  SUNMemory Brecv_;
  SUNMemory Frecv_;

  SUNMemoryType memoryType()
  {
    SUNMemory test;
    if (!SUNMemoryHelper_Alloc(memhelp, &test, sizeof(REAL), SUNMEMTYPE_PINNED,
                               nullptr))
    {
      SUNMemoryHelper_Dealloc(memhelp, test, nullptr);
      return (SUNMEMTYPE_PINNED);
    }
    if (!SUNMemoryHelper_Alloc(memhelp, &test, sizeof(REAL), SUNMEMTYPE_DEVICE,
                               nullptr))
    {
      SUNMemoryHelper_Dealloc(memhelp, test, nullptr);
      return (SUNMEMTYPE_DEVICE);
    }
    if (!SUNMemoryHelper_Alloc(memhelp, &test, sizeof(REAL), SUNMEMTYPE_UVM,
                               nullptr))
    {
      SUNMemoryHelper_Dealloc(memhelp, test, nullptr);
      return (SUNMEMTYPE_UVM);
    }
    else { return (SUNMEMTYPE_HOST); }
  }
};

} // namespace sundials_tools

#endif
