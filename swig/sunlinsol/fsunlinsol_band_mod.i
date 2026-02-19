// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
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

%module fsunlinsol_band_mod

// include code common to all nvector implementations
%include "fsunlinsol.i"

%{
#include "sunlinsol/sunlinsol_band.h"
%}

// sunlinsol_impl macro defines some ignore and inserts with the linear solver name appended
%sunlinsol_impl(Band)

// Process and wrap functions in the following files
%include "sunlinsol/sunlinsol_band.h"

