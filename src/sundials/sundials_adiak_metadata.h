/* -----------------------------------------------------------------
 * Programmer(s): Yu Pan and Cody Balos @ LLNL
 * -----------------------------------------------------------------
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
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_ADIAK_METADATA_H
#define _SUNDIALS_ADIAK_METADATA_H

#ifdef SUNDIALS_ADIAK_ENABLED
#include <adiak.h>

static void sunAdiakCollectMetadata()
{
  adiak_launchdate();
  adiak_executable();
  adiak_cmdline();
  adiak_clustername();

  adiak_job_size();
  adiak_num_hosts();

  adiak_namevalue("c_compiler", 2, NULL, "%s", SUN_C_COMPILER);
  adiak_namevalue("c_compiler_version", 2, NULL, "%s", SUN_C_COMPILER_VERSION);
  adiak_namevalue("c_compiler_flags", 2, NULL, "%s", SUN_C_COMPILER_FLAGS);

  adiak_namevalue("cxx_compiler", 2, NULL, "%s", SUN_CXX_COMPILER);
  adiak_namevalue("cxx_compiler_version", 2, NULL, "%s",
                  SUN_CXX_COMPILER_VERSION);
  adiak_namevalue("cxx_compiler_flags", 2, NULL, "%s", SUN_CXX_COMPILER_FLAGS);

  adiak_namevalue("fortran_compiler", 2, NULL, "%s", SUN_FORTRAN_COMPILER);
  adiak_namevalue("fortran_compiler_version", 2, NULL, "%s",
                  SUN_FORTRAN_COMPILER_VERSION);
  adiak_namevalue("fortran_compiler_flags", 2, NULL, "%s",
                  SUN_FORTRAN_COMPILER_FLAGS);

  adiak_namevalue("sundials_version", 2, NULL, "%s", SUNDIALS_VERSION);
  adiak_namevalue("sundials_git_version", 2, NULL, "%s", SUNDIALS_GIT_VERSION);
  adiak_namevalue("build_type", 2, NULL, "%s", SUN_BUILD_TYPE);
  adiak_namevalue("third_party_libraries", 2, NULL, "%s", SUN_TPL_LIST);
#ifdef SUN_JOB_ID
  adiak_namevalue("job_id", 2, NULL, "%s", SUN_JOB_ID);
#endif
  adiak_namevalue("job_start_time", 2, NULL, "%s", SUN_JOB_START_TIME);

#ifdef SUNDIALS_SPACK_VERSION
  adiak_namevalue("spack_version", 2, NULL, "%s", SUNDIALS_SPACK_VERSION);
#endif

#ifdef SUNDIALS_GINKGO_ENABLED
  adiak_namevalue("ginkgo_version", 2, NULL, "%s", SUN_GINKGO_VERSION);
#endif

#ifdef SUNDIALS_HYPRE_ENABLED
  adiak_namevalue("hypre_version", 2, NULL, "%s", SUN_HYPRE_VERSION);
#endif

#ifdef SUNDIALS_KLU_ENABLED
  adiak_namevalue("klu_version", 2, NULL, "%s", SUN_KLU_VERSION);
#endif

#ifdef SUNDIALS_KOKKOS_ENABLED
  adiak_namevalue("kokkos_version", 2, NULL, "%s", SUN_KOKKOS_VERSION);
#endif

#ifdef SUNDIALS_KOKKOS_KERNELS_ENABLED
  adiak_namevalue("kokkos_kernels_version", 2, NULL, "%s",
                  SUN_KOKKOS_KERNELS_VERSION);
#endif

#ifdef SUNDIALS_BLAS_LAPACK_ENABLED
  adiak_namevalue("lapack_version", 2, NULL, "%s", SUN_LAPACK_VERSION);
#endif

#ifdef SUNDIALS_MAGMA_ENABLED
  adiak_namevalue("magma_version", 2, NULL, "%s", SUN_MAGMA_VERSION);
#endif

#if SUNDIALS_MPI_ENABLED
  adiak_namevalue("mpi_c_compiler", 2, NULL, "%s", SUN_MPI_C_COMPILER);
  adiak_namevalue("mpi_c_version", 2, NULL, "%s", SUN_MPI_C_VERSION);

  adiak_namevalue("mpi_cxx_compiler", 2, NULL, "%s", SUN_MPI_CXX_COMPILER);
  adiak_namevalue("mpi_cxx_version", 2, NULL, "%s", SUN_MPI_CXX_VERSION);

  adiak_namevalue("mpi_fortran_compiler", 2, NULL, "%s",
                  SUN_MPI_FORTRAN_COMPILER);
  adiak_namevalue("mpi_fortran_version", 2, NULL, "%s", SUN_MPI_FORTRAN_VERSION);
#endif

#ifdef SUNDIALS_ONEMKL_ENABLED
  adiak_namevalue("onemkl_version", 2, NULL, "%s", SUN_ONEMKL_VERSION);
#endif

#ifdef SUNDIALS_OPENMP_ENABLED
  adiak_namevalue("openmp_version", 2, NULL, "%s", SUN_OPENMP_VERSION);
#endif

#ifdef SUNDIALS_PETSC_ENABLED
  adiak_namevalue("petsc_version", 2, NULL, "%s", SUN_PETSC_VERSION);
#endif

#ifdef SUNDIALS_PTHREADS_ENABLED
  adiak_namevalue("pthreads_version", 2, NULL, "%s", SUN_PTHREADS_VERSION);
#endif

#ifdef SUNDIALS_RAJA_ENABLED
  adiak_namevalue("raja_version", 2, NULL, "%s", SUN_RAJA_VERSION);
#endif

#ifdef SUNDIALS_SUPERLUDIST_ENABLED
  adiak_namevalue("superludist_version", 2, NULL, "%s", SUN_SUPERLUDIST_VERSION);
#endif

#ifdef SUNDIALS_SUPERLUMT_ENABLED
  adiak_namevalue("superlumt_version", 2, NULL, "%s", SUN_SUPERLUMT_VERSION);
#endif

#ifdef SUNDIALS_TRILLINOS_ENABLED
  adiak_namevalue("trillinos_version", 2, NULL, "%s", SUN_TRILLINOS_VERSION);
#endif

#ifdef SUNDIALS_XBRAID_ENABLED
  adiak_namevalue("xbraid_version", 2, NULL, "%s", SUN_XBRAID_VERSION);
#endif

#ifdef SUNDIALS_CUDA_ENABLED
  adiak_namevalue("cuda_version", 2, NULL, "%s", SUN_CUDA_VERSION);
  adiak_namevalue("cuda_compiler", 2, NULL, "%s", SUN_CUDA_COMPILER);
  adiak_namevalue("cuda_architectures", 2, NULL, "%s", SUN_CUDA_ARCHITECTURES);
#endif

#ifdef SUNDIALS_HIP_ENABLED
  adiak_namevalue("hip_version", 2, NULL, "%s", SUN_HIP_VERSION);
  adiak_namevalue("amdgpu_targets", 2, NULL, "%s", SUN_AMDGPU_TARGETS);
#endif
}
#endif

#endif
