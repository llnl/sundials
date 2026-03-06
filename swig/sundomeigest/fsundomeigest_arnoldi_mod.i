// ---------------------------------------------------------------
// Programmer: Mustafa Aggul @ SMU
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

%module fsundomeigest_arnoldi_mod

// include code common to all implementations
%include "fsundomeigest.i"

// Ignore command-line processing functions since they are not supported in Fortran
%ignore SUNDomEigEstimator_SetOptions_Arnoldi;

%{
#include "sundomeigest/sundomeigest_arnoldi.h"
%}

// sundomeigest_impl macro defines some ignore and inserts with the estimator name appended
%sundomeigest_impl(ARNOLDI)

// Process and wrap functions in the following files
%include "sundomeigest/sundomeigest_arnoldi.h"