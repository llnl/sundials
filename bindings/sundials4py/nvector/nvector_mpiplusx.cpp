/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <mpi.h>
#include <nvector/nvector_mpiplusx.h>
#include <sundials/sundials_nvector.hpp>

#include "nvector_array_helpers.hpp"

namespace nb = nanobind;

namespace sundials4py {

void bind_nvector_mpiplusx(nb::module_& m)
{
  m.def(
    "N_VMake_MPIPlusX",
    [](int fortran_comm, N_Vector local, SUNContext sunctx)
    {
      MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(fortran_comm));
      return nvector_detail::wrap_nvector(N_VMake_MPIPlusX(comm, local, sunctx));
    },
    nb::arg("fortran_comm"), nb::arg("local"), nb::arg("sunctx"),
    nb::keep_alive<0, 2>(), nb::keep_alive<0, 3>());

  m.def(
    "N_VGetLocalVector_MPIPlusX",
    [](N_Vector v)
    {
      auto* local = N_VGetLocalVector_MPIPlusX(v);
      return std::shared_ptr<std::remove_pointer_t<N_Vector>>(local,
                                                              [](auto*) {});
    },
    nb::arg("v"), nb::keep_alive<0, 1>());
}

} // namespace sundials4py
