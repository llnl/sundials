// ---------------------------------------------------------------
// Programmer: Cody J. Balos @ LLNL
// ---------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2025, Lawrence Livermore National Security,
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

%module farkode_mristep_mod

%include "../sundials/fsundials.i"

// include the header file(s) in the c wrapper that is generated
%{
#include "arkode/arkode_mristep.h"
#include "arkode/arkode_mristep_deprecated.h"
%}

// Load the typedefs and generate a "use" statements in the module
%import "farkode_mod.i"

// Treat MRIStepInnerStepper and MRIStepCoupling as void*
%apply void* { MRIStepInnerStepper };
%apply void* { MRIStepCoupling };

// Process definitions from these files
%include "arkode/arkode_mristep.h"
%include "arkode/arkode_mristep_deprecated.h"

