# ----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ----------------------------------------------------------------
doc_version = "v7.4.0"
sundials_version = "v7.4.0"
arkode_version = "v6.4.0"
cvode_version = "v7.4.0"
cvodes_version = "v7.4.0"
ida_version = "v7.4.0"
idas_version = "v6.4.0"
kinsol_version = "v7.4.0"
year = "2025"

# Warn about all references where the target cannot be found
nitpicky = True

# List of tuples (type, target) to ignore when generating "nitpicky" warnings
nitpick_ignore = [
    # C/C++ types Sphinx does not seem to know about
    ("c:identifier", "FILE"),
    ("cpp:identifier", "FILE"),
    ("c:identifier", "size_t"),
    ("cpp:identifier", "size_t"),
    ("c:identifier", "int64_t"),
    # CUDA
    ("cpp:identifier", "cudaStream_t"),
    ("c:identifier", "cusparseHandle_t"),
    ("c:identifier", "cusparseMatDescr_t"),
    ("c:identifier", "cusolverHandle_t"),
    # Ginkgo
    ("cpp:identifier", "gko"),
    ("cpp:identifier", "gko::dim<2>"),
    ("cpp:identifier", "gko::Executor"),
    ("cpp:identifier", "gko::LinOp"),
    # HIP
    ("cpp:identifier", "hipStream_t"),
    # hypre
    ("c:identifier", "hypre_ParVector"),
    # Kokkos
    ("cpp:identifier", "ExecutionSpace"),
    ("cpp:identifier", "ExecutionSpace::memory_space"),
    ("cpp:identifier", "Kokkos"),
    ("cpp:identifier", "Kokkos::DefaultExecutionSpace"),
    ("cpp:identifier", "Kokkos::RangePolicy<exec_space>"),
    ("cpp:identifier", "Kokkos::View<sunrealtype***, memory_space>"),
    ("cpp:identifier", "Kokkos::View<sunrealtype*, MemorySpace>"),
    ("cpp:identifier", "Kokkos::Rank<3>"),
    ("cpp:identifier", "Kokkos::TeamPolicy<exec_space>"),
    ("cpp:identifier", "Kokkos::TeamPolicy<exec_space>::member_type"),
    ("cpp:identifier", "Kokkos::MDRangePolicy<exec_space, Kokkos::Rank<3>>"),
    # MPI
    ("c:identifier", "MPI_Comm"),
    # PETSc
    ("c:identifier", "SNES"),
    ("c:identifier", "PetscErrorCode"),
    ("c:identifier", "Vec"),
    # SuperLU
    ("c:identifier", "gridinfo_t"),
    ("c:identifier", "SuperMatrix"),
    ("c:identifier", "gridinfo_t"),
    ("c:identifier", "xLUstruct_t"),
    ("c:identifier", "xScalePermstruct_t"),
    ("c:identifier", "xSOLVEstruct_t"),
    ("c:identifier", "SuperLUStat_t"),
    ("c:identifier", "superlu_dist_options_t"),
    # SYCL
    ("cpp:identifier", "sycl"),
    ("cpp:identifier", "sycl::queue"),
    # Trilinos
    ("cpp:identifier", "Tpetra"),
    ("cpp:identifier", "Tpetra::Vector<sunrealtype, int, sunindextype>"),
    ("cpp:identifier", "Teuchos"),
    ("cpp:identifier", "Teuchos::RCP<vector_type>"),
    # XBraid
    ("c:identifier", "braid_AccessStatus"),
    ("c:identifier", "braid_App"),
    ("c:identifier", "braid_BufferStatus"),
    ("c:identifier", "braid_Core"),
    ("c:identifier", "braid_Int"),
    ("c:identifier", "braid_PtFcnAccess"),
    ("c:identifier", "braid_PtFcnInit"),
    ("c:identifier", "braid_PtFcnSpatialNorm"),
    ("c:identifier", "braid_PtFcnStep"),
    ("c:identifier", "braid_Real"),
    ("c:identifier", "braid_StepStatus"),
    ("c:identifier", "braid_Vector"),
    # C types referenced in C++ functions, not sure how to fix
    ("cpp:identifier", "sunbooleantype"),
    ("cpp:identifier", "suncountertype"),
    ("cpp:identifier", "sunindextype"),
    ("cpp:identifier", "sunrealtype"),
    ("cpp:identifier", "SUNErrCode"),
    ("cpp:identifier", "SUNContext"),
    ("cpp:identifier", "N_Vector"),
    ("cpp:identifier", "SUNMatrix"),
    ("cpp:identifier", "SUNLinearSolver"),
    ("cpp:identifier", "SUNMemoryHelper"),
    ("cpp:identifier", "SUNMemoryType"),
    # C++ namespaces don't seem to work as expected, not sure how to fix
    ("cpp:identifier", "sundials"),
    ("cpp:identifier", "sundials::cuda"),
    ("cpp:identifier", "sundials::hip"),
    ("cpp:identifier", "sundials::sycl"),
    ("cpp:identifier", "sundials::ginkgo"),
    # Experimental or internal namespaces and types
    ("cpp:identifier", "sundials::impl"),
    ("cpp:identifier", "sundials::impl::BaseNVector"),
    ("cpp:identifier", "sundials::impl::BaseMatrix"),
    ("cpp:identifier", "sundials::impl::BaseLinearSolver"),
    ("cpp:identifier", "sundials::ConvertibleTo<N_Vector>"),
    ("cpp:identifier", "sundials::ConvertibleTo<SUNMatrix>"),
    ("cpp:identifier", "sundials::ConvertibleTo<SUNLinearSolver>"),
    # Defined types in Kokkos vector that are not found by Sphinx
    ("cpp:identifier", "view_type::size_type"),
    ("cpp:identifier", "view_type::HostMirror"),
    # Template parameter that causes an error in Kokkos matrix
    ("cpp:identifier", "MatrixType"),
    # C++ types referenced in "C" functions, should probably switch
    # documentation to use .. cpp:function rather than .. c:function
    ("c:identifier", "SUNCudaExecPolicy"),
    ("c:identifier", "SUNHipExecPolicy"),
]
