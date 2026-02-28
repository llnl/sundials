// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2025-2026, Lawrence Livermore National Security,
// University of Maryland Baltimore County, and the SUNDIALS contributors.
// Copyright (c) 2013-2025, Lawrence Livermore National Security
// and Southern Methodist University.
// Copyright (c) 2002-2013, Lawrence Livermore National Security.
// All rights reserved.
//
// See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-3-Clause
// SUNDIALS Copyright End
// ---------------------------------------------------------------
// Swig interface file
// ---------------------------------------------------------------

%module fsundials_core_mod

%include "fsundials_types.i"
%include "fsundials.i"
%include "fsundials_context.i"
%include "fsundials_profiler.i"
%include "fsundials_logger.i"
%include "fsundials_futils.i"
%include "fsundials_nvector.i"
%include "fsundials_matrix.i"
%include "fsundials_linearsolver.i"
%include "fsundials_nonlinearsolver.i"
%include "fsundials_adaptcontroller.i"
%include "fsundials_stepper.i"
%include "fsundials_memory.i"
%include "fsundials_adjoint.i"
%include "fsundials_domeigestimator.i"
%include "fcopyright.i"
