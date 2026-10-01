/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ UMBC
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * Umbrella header included by every sundials4py translation unit.
 *
 * It supplies the basic types and exceptions, the four Python-subclassable
 * "custom" object base classes, and the nanobind type casters that make those
 * subclasses acceptable wherever a raw SUNDIALS pointer is expected. Because the
 * casters must be visible in every translation unit that binds a SUNDIALS
 * signature -- otherwise different units would disagree about how the type is
 * passed -- they are reached through this one header rather than included ad hoc.
 *
 * The implementations live in the per-family headers listed below.
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_TYPES_HPP
#define _SUNDIALS4PY_TYPES_HPP

/* Basic array typedefs and the sundials4py exception hierarchy. */
#include "sundials4py_core_types.hpp"

/* Machinery shared by the custom object families (callback revocation, shutdown
   safety, override detection, tagged content). */
#include "sundials4py_custom_object.hpp"

/* The custom object base classes themselves. */
#include "sundials_adaptcontroller_custom.hpp"
#include "sundials_linearsolver_custom.hpp"
#include "sundials_matrix_custom.hpp"
#include "sundials_nonlinearsolver_custom.hpp"

/* type_caster specializations for _generic_SUNMatrix and friends. Must come last,
   since it names all four classes above. */
#include "sundials4py_custom_casters.hpp"

#endif // _SUNDIALS4PY_TYPES_HPP
