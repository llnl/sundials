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
 * SUNDIALS N_Vector MPIPlusX bindings. Hand-written code for functions
 * that need special treatment, plus the generated code produced with the
 * generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <nvector/nvector_mpiplusx.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_nvector.hpp>

#include "../sundials/sundials_comm.hpp"
#include "nvector_array_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_mpiplusx(nb::module_& m)
{
#include "nvector_mpiplusx_generated.hpp"

  m.def(
    "N_VMake_MPIPlusX",
    [](nb::handle comm, N_Vector local,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      N_Vector v = N_VMake_MPIPlusX(to_suncomm(comm), local, sunctx);
      if (!v) { throw sundials4py::error_returned("N_VMake_MPIPlusX failed"); }
      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(v);
    },
    nb::arg("comm").none(), nb::arg("local"), nb::arg("sunctx"),
    nb::keep_alive<0, 2>() /* keep the local vector alive as long as the N_Vector is */,
    nb::keep_alive<0, 3>() /* keep the SUNContext alive as long as the N_Vector is */,
    "Create an MPIPlusX N_Vector wrapping a local vector.\n\n"
    "This call is collective over comm (it performs an MPI_Allreduce). comm "
    "may be None or SUN_COMM_NULL, in which case no communication occurs. "
    "The vector duplicates comm internally.");

  m.def(
    "N_VGetLocalVector_MPIPlusX",
    [](nb::handle v_obj) -> nb::object
    {
      N_Vector v;
      if (!nb::try_cast<N_Vector>(v_obj, v))
      {
        throw nb::type_error("v must be an N_Vector");
      }
      if (N_VGetVectorID(v) != SUNDIALS_NVEC_MPIPLUSX)
      {
        throw sundials4py::illegal_value("v is not an MPIPlusX N_Vector");
      }
      N_Vector local = N_VGetLocalVector_MPIPlusX(v);
      if (!local)
      {
        throw sundials4py::error_returned(
          "N_VGetLocalVector_MPIPlusX returned NULL");
      }
      return nvector_detail::child_vector_object(v_obj, local);
    },
    nb::arg("v"),
    "Return the local vector of an MPIPlusX N_Vector.\n\n"
    "If the local vector already has a Python wrapper (for example the one "
    "passed to N_VMake_MPIPlusX), that original wrapper is returned. "
    "Otherwise a non-owning reference wrapper that keeps v alive is returned. "
    "Do not keep it beyond the lifetime of v (for example, do not store it "
    "from inside a callback).");

  m.def(
    "N_VGetLocalLength_MPIPlusX",
    [](N_Vector v) { return N_VGetLocalLength_MPIPlusX(v); }, nb::arg("v"));
}

} // namespace sundials4py
