# ----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# ----------------------------------------------------------------
doc_version = "develop"
sundials_version = "v7.5.0"
arkode_version = "v6.5.0"
cvode_version = "v7.5.0"
cvodes_version = "v7.5.0"
ida_version = "v7.5.0"
idas_version = "v6.5.0"
kinsol_version = "v7.5.0"
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
    ("cpp:identifier", "gko::size_type"),
    ("cpp:identifier", "gko::batch_dim<2>"),
    ("cpp:identifier", "gko::batch"),
    ("cpp:identifier", "gko::batch::stop"),
    ("cpp:identifier", "gko::batch::stop::tolerance_type"),
    ("cpp:identifier", "gko::batch::BatchLinOpFactory"),
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
    # Python
    ("py:class", "typing_extensions.CapsuleType"),
    ("py:class", "types.CapsuleType"),
    ("py:class", "collections.abc.Callable"),
    ("py:class", "collections.abc.Sequence"),
    ("py:class", "numpy.ndarray[dtype=float64, shape=(*), order='C']"),
    ("py:class", "numpy.ndarray[dtype=float64, shape=(*, *), order='F']"),
    ("py:class", "numpy.ndarray[dtype=int64, shape=(*), order='C']"),
    # C functions we don't document explicitly
    ("c:func", "ARKodePrintMem"),
    ("c:func", "ARKodeSetPostprocessStageFn"),
    ("c:func", "ARKodeSetPostprocessStepFn"),
    ("c:func", "CVodeGetAdjDataPointHermite"),
    ("c:func", "CVodeGetAdjDataPointPolynomial"),
    ("c:func", "CVodeQuadSStolerancesB"),
    ("c:func", "CVodeQuadSVtolerancesB"),
    ("c:func", "CVodeSetConstraintsB"),
    ("c:func", "CVodeSetInitStepB"),
    ("c:func", "CVodeSetMaxNumStepsB"),
    ("c:func", "CVodeSetMaxOrdB"),
    ("c:func", "CVodeSetMaxStepB"),
    ("c:func", "CVodeSetMinStepB"),
    ("c:func", "CVodeSetQuadErrConB"),
    ("c:func", "CVodeSetStabLimDetB"),
    ("c:func", "IDAGetAdjDataPointHermite"),
    ("c:func", "IDAGetAdjDataPointPolynomial"),
    ("c:func", "IDAQuadSStolerancesB"),
    ("c:func", "IDAQuadSVtolerancesB"),
    ("c:func", "IDASetConstraintsB"),
    ("c:func", "IDASetIdB"),
    ("c:func", "IDASetInitStepB"),
    ("c:func", "IDASetMaxNumStepsB"),
    ("c:func", "IDASetMaxOrdB"),
    ("c:func", "IDASetMaxStepB"),
    ("c:func", "IDASetQuadErrConB"),
    ("c:func", "IDASetSuppressAlgB"),
    ("c:func", "N_VEnableDotProdMultiLocal_ManyVector"),
    ("c:func", "SUNAdaptController_GetFastController_MRIHTol"),
    ("c:func", "SUNAdaptController_GetSlowController_MRIHTol"),
    ("c:func", "SUNClassicalGS"),
    ("c:func", "SUNMemoryHelper_ImplementsRequiredOps"),
    ("c:func", "SUNModifiedGS"),
    ("c:func", "SUNNonlinSol_FixedPointSens"),
    ("c:func", "SUNNonlinSol_NewtonSens"),
    ("c:func", "SUNQRAdd_CGS2"),
    ("c:func", "SUNQRAdd_DCGS2_SB"),
    ("c:func", "SUNQRAdd_DCGS2"),
    ("c:func", "SUNQRAdd_ICWY_SB"),
    ("c:func", "SUNQRAdd_ICWY"),
    ("c:func", "SUNQRAdd_MGS"),
    ("c:func", "SUNQRfact"),
    ("c:func", "SUNQRsol"),
    ("c:func", "SUNSparseMatrix_ToCSC"),
    ("c:func", "SUNSparseMatrix_ToCSR"),
]
