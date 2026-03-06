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

// Include shared configuration
%include "../sundials/fsundials.i"

%{
#include "sundials/sundials_domeigestimator.h"
%}

// Load the typedefs and generate "use" statements
%import "../sundials/fsundials_core_mod.i"

// Macro for creating an interface to an SUNDomEigEstimator
%define %sundomeigest_impl(TYPE)
  %ignore SUNDomEigEstimatorContent_## TYPE ##;
%enddef
