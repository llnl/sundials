# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Setup third-party libraries
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Setup MPI, OpenMP, and OpenMP offload first as other TPLs may need targets or
# variables corresponding to these TPLs.
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Find MPI
# ---------------------------------------------------------------

if(ENABLE_MPI)
  include(SundialsMPI)
  list(APPEND SUNDIALS_TPL_LIST "MPI")
endif()

# ---------------------------------------------------------------
# Find OpenMP
# ---------------------------------------------------------------

if(ENABLE_OPENMP)
  include(SundialsOpenMP)
  list(APPEND SUNDIALS_TPL_LIST "OPENMP")
endif()

# ---------------------------------------------------------------
# Find OpenMP with device offloading
# --------------------------------------------------------------

if(ENABLE_OPENMP_DEVICE)
  include(SundialsOpenMP)
  list(APPEND SUNDIALS_TPL_LIST "OPENMP_DEVICE")
endif()

# ---------------------------------------------------------------
# Setup other TPLs (listed in alphabetical order)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Find (and test) the Adiak libraries
# ---------------------------------------------------------------

if(ENABLE_ADIAK)
  include(SundialsAdiak)
  list(APPEND SUNDIALS_TPL_LIST "ADIAK")
endif()

# ---------------------------------------------------------------
# Find (and test) the Caliper libraries
# ---------------------------------------------------------------

if(ENABLE_CALIPER)
  include(SundialsCaliper)
  list(APPEND SUNDIALS_TPL_LIST "CALIPER")
endif()

# ---------------------------------------------------------------
# Find (and test) the Ginkgo libraries
# ---------------------------------------------------------------

if(ENABLE_GINKGO)
  include(SundialsGinkgo)
  list(APPEND SUNDIALS_TPL_LIST "GINKGO")
endif()

# ---------------------------------------------------------------
# Find (and test) the hypre libraries
# ---------------------------------------------------------------

if(ENABLE_HYPRE)
  include(SundialsHypre)
  list(APPEND SUNDIALS_TPL_LIST "HYPRE")
endif()

# ---------------------------------------------------------------
# Find (and test) Kokkos
# ---------------------------------------------------------------

if(ENABLE_KOKKOS)
  include(SundialsKokkos)
  list(APPEND SUNDIALS_TPL_LIST "KOKKOS")
endif()

# ---------------------------------------------------------------
# Find (and test) Kokkos Kernels
# ---------------------------------------------------------------

if(ENABLE_KOKKOS_KERNELS)
  include(SundialsKokkosKernels)
  list(APPEND SUNDIALS_TPL_LIST "KOKKOS_KERNELS")
endif()

# ---------------------------------------------------------------
# Find (and test) the KLU libraries
# ---------------------------------------------------------------

if(ENABLE_KLU)
  include(SundialsKLU)
  list(APPEND SUNDIALS_TPL_LIST "KLU")
endif()

# ---------------------------------------------------------------
# Find (and test) the LAPACK and BLAS libraries
# ---------------------------------------------------------------

if(ENABLE_LAPACK)
  include(SundialsLapack)
  list(APPEND SUNDIALS_TPL_LIST "BLAS_LAPACK")
endif()

# ---------------------------------------------------------------
# Find (and test) the MAGMA libraries
# ---------------------------------------------------------------

if(ENABLE_MAGMA)
  include(SundialsMAGMA)
  list(APPEND SUNDIALS_TPL_LIST "MAGMA")
endif()

# ---------------------------------------------------------------
# Find (and test) the oneMKL libraries
# ---------------------------------------------------------------

if(ENABLE_ONEMKL)
  include(SundialsONEMKL)
  list(APPEND SUNDIALS_TPL_LIST "ONEMKL")
endif()

# ---------------------------------------------------------------
# Find (and test) the PETSc libraries
# ---------------------------------------------------------------

if(ENABLE_PETSC)
  include(SundialsPETSC)
  list(APPEND SUNDIALS_TPL_LIST "PETSC")
endif()

# ---------------------------------------------------------------
# Find PThreads
# ---------------------------------------------------------------

if(ENABLE_PTHREAD)
  include(SundialsPthread)
  list(APPEND SUNDIALS_TPL_LIST "PTHREAD")
endif()

# -------------------------------------------------------------
# Find (and test) RAJA
# -------------------------------------------------------------

if(ENABLE_RAJA)
  include(SundialsRAJA)
  list(APPEND SUNDIALS_TPL_LIST "RAJA")
endif()

# ---------------------------------------------------------------
# Find (and test) the SuperLUDIST libraries
# ---------------------------------------------------------------

if(ENABLE_SUPERLUDIST)
  include(SundialsSuperLUDIST)
  list(APPEND SUNDIALS_TPL_LIST "SUPERLUDIST")
endif()

# ---------------------------------------------------------------
# Find (and test) the SUPERLUMT libraries
# ---------------------------------------------------------------

if(ENABLE_SUPERLUMT)
  include(SundialsSuperLUMT)
  list(APPEND SUNDIALS_TPL_LIST "SUPERLUMT")
endif()

# -------------------------------------------------------------
# Find (and test) Trilinos
# -------------------------------------------------------------

if(ENABLE_TRILINOS)
  include(SundialsTrilinos)
  list(APPEND SUNDIALS_TPL_LIST "TRILINOS")
endif()

# -------------------------------------------------------------
# Find (and test) XBraid
# -------------------------------------------------------------

if(ENABLE_XBRAID)
  include(SundialsXBRAID)
  list(APPEND SUNDIALS_TPL_LIST "XBRAID")
endif()
