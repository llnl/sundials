/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------
 * Conversion between Python communicator objects and SUNComm.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_SUNDIALS4PY_COMM_HPP
#define SUNDIALS_SUNDIALS4PY_COMM_HPP

#include <sundials/sundials_types.h>

#include "sundials4py.hpp"

namespace sundials4py {

/// Converts a Python object to a SUNComm. Accepts None or the int 0 (the
/// null communicator) and, in MPI builds, an mpi4py.MPI.Comm. Anything else
/// raises TypeError. The returned handle is borrowed from the Python object.
SUNComm to_suncomm(nanobind::handle obj);

/// Converts a SUNComm to a Python object. The null communicator becomes the
/// int 0. In MPI builds, any other communicator becomes a non-owning
/// mpi4py.MPI.Comm view that must not be freed by the caller.
nanobind::object from_suncomm(SUNComm c);

} // namespace sundials4py

#endif // SUNDIALS_SUNDIALS4PY_COMM_HPP
