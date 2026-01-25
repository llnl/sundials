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

%module fcvode_mod

%include "../sundials/fsundials.i"

// Ignore command-line processing functions since they are not supported in Fortran
%ignore CVodeSetOptions;

%{
#include "cvode/cvode.h"
#include "cvode/cvode_bandpre.h"
#include "cvode/cvode_bbdpre.h"
#include "cvode/cvode_diag.h"
#include "cvode/cvode_ls.h"
#include "cvode/cvode_proj.h"
%}

%import "../sundials/fsundials_core_mod.i"

// Process definitions from these files
%include "cvode/cvode.h"
%include "cvode/cvode_bandpre.h"
%include "cvode/cvode_bbdpre.h"
%include "cvode/cvode_diag.h"
%include "cvode/cvode_ls.h"
%include "cvode/cvode_proj.h"
