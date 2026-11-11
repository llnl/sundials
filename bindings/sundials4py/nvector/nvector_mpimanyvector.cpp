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
 * SUNDIALS N_Vector MPIManyVector bindings. Hand-written code for
 * functions that need special treatment, plus the generated code produced
 * with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <nvector/nvector_mpimanyvector.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_nvector.hpp>

#include "../sundials/sundials_comm.hpp"
#include "nvector_array_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

// Validates the subvector sequence and returns it as (a) a tuple that the new
// vector keeps alive (not the caller's mutable list) and (b) raw pointers.
nb::tuple validate_subvectors(sunindextype num_subvectors, nb::sequence subvecs,
                              std::vector<N_Vector>& out)
{
  if (num_subvectors < 1)
  {
    throw sundials4py::illegal_value("num_subvectors must be at least 1");
  }
  if (static_cast<sunindextype>(nb::len(subvecs)) != num_subvectors)
  {
    throw sundials4py::illegal_value(
      "len(subvecs) must be equal to num_subvectors");
  }

  nb::list items;
  for (sunindextype i = 0; i < num_subvectors; i++)
  {
    nb::object item = subvecs[nb::int_(i)];
    if (item.is_none())
    {
      throw sundials4py::illegal_value("subvecs must not contain None");
    }
    N_Vector sv;
    if (!nb::try_cast<N_Vector>(item, sv))
    {
      throw nb::type_error("subvecs must contain only N_Vectors");
    }
    out.push_back(sv);
    items.append(item);
  }
  return nb::tuple(items);
}

nb::object wrap_new_mpimanyvector(N_Vector v, nb::tuple keep, nb::handle sunctx)
{
  nb::object result = nb::cast(
    our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(v));
  nvector_detail::keep_alive_objects(result.ptr(), keep.ptr());
  nvector_detail::keep_alive_objects(result.ptr(), sunctx.ptr());
  return result;
}

} // namespace

void bind_nvector_mpimanyvector(nb::module_& m)
{
#include "nvector_mpimanyvector_generated.hpp"

  m.def(
    "N_VMake_MPIManyVector",
    [](nb::handle comm, sunindextype num_subvectors, nb::sequence subvecs,
       nb::handle sunctx_obj) -> nb::object
    {
      SUNContext sunctx;
      if (!nb::try_cast<SUNContext>(sunctx_obj, sunctx))
      {
        throw nb::type_error("sunctx must be a SUNContext");
      }
      std::vector<N_Vector> vecs;
      nb::tuple keep = validate_subvectors(num_subvectors, subvecs, vecs);
      N_Vector v     = N_VMake_MPIManyVector(to_suncomm(comm), num_subvectors,
                                             vecs.data(), sunctx);
      if (!v)
      {
        throw sundials4py::error_returned("N_VMake_MPIManyVector failed");
      }
      return wrap_new_mpimanyvector(v, keep, sunctx_obj);
    },
    nb::arg("comm").none(), nb::arg("num_subvectors"), nb::arg("subvecs"),
    nb::arg("sunctx"),
    "Create an MPIManyVector from subvectors.\n\n"
    "This call is collective over comm (it performs an MPI_Allreduce). comm "
    "may be None or SUN_COMM_NULL, in which case no communication occurs. "
    "The vector keeps the subvectors and the SUNContext alive.");

  m.def(
    "N_VNew_MPIManyVector",
    [](sunindextype num_subvectors, nb::sequence subvecs,
       nb::handle sunctx_obj) -> nb::object
    {
      SUNContext sunctx;
      if (!nb::try_cast<SUNContext>(sunctx_obj, sunctx))
      {
        throw nb::type_error("sunctx must be a SUNContext");
      }
      std::vector<N_Vector> vecs;
      nb::tuple keep = validate_subvectors(num_subvectors, subvecs, vecs);
      N_Vector v = N_VNew_MPIManyVector(num_subvectors, vecs.data(), sunctx);
      if (!v)
      {
        throw sundials4py::error_returned(
          "N_VNew_MPIManyVector failed; at least one subvector must have a "
          "communicator (for example an MPIPlusX vector), otherwise use "
          "N_VMake_MPIManyVector");
      }
      return wrap_new_mpimanyvector(v, keep, sunctx_obj);
    },
    nb::arg("num_subvectors"), nb::arg("subvecs"), nb::arg("sunctx"),
    "Create an MPIManyVector whose communicator is taken from the "
    "subvectors, at least one of which must be MPI-aware.\n\n"
    "This call is collective. The vector keeps the subvectors and the "
    "SUNContext alive.");

  m.def(
    "N_VGetSubvector_MPIManyVector",
    [](nb::handle v_obj, sunindextype vec_num) -> nb::object
    {
      N_Vector v;
      if (!nb::try_cast<N_Vector>(v_obj, v))
      {
        throw nb::type_error("v must be an N_Vector");
      }
      if (N_VGetVectorID(v) != SUNDIALS_NVEC_MPIMANYVECTOR)
      {
        throw sundials4py::illegal_value("v is not an MPIManyVector N_Vector");
      }
      if (vec_num < 0 || vec_num >= N_VGetNumSubvectors_MPIManyVector(v))
      {
        throw nb::index_error("subvector index out of range");
      }
      N_Vector sub = N_VGetSubvector_MPIManyVector(v, vec_num);
      if (!sub)
      {
        throw sundials4py::error_returned(
          "N_VGetSubvector_MPIManyVector returned NULL");
      }
      return nvector_detail::child_vector_object(v_obj, sub);
    },
    nb::arg("v"), nb::arg("vec_num"),
    "Return subvector vec_num of an MPIManyVector.\n\n"
    "If the subvector already has a Python wrapper (for example one passed "
    "at construction), that original wrapper is returned. Otherwise a "
    "non-owning reference wrapper that keeps v alive is returned. Do not "
    "keep it beyond the lifetime of v (for example, do not store it from "
    "inside a callback).");
}

} // namespace sundials4py
