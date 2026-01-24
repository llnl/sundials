/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
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
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_core.hpp>
#include <sunmatrix/sunmatrix_dense.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunmatrix_dense(nb::module_& m)
{
#include "sunmatrix_dense_generated.hpp"

  m.def(
    "SUNDenseMatrix_Data",
    [](SUNMatrix A)
    {
      auto rows  = static_cast<size_t>(SUNDenseMatrix_Rows(A));
      auto cols  = static_cast<size_t>(SUNDenseMatrix_Columns(A));
      auto owner = nb::find(A);
      auto ptr   = SUNDenseMatrix_Data(A);
      // SUNDenseMatrix_Data returns a column-major ordered array (i.e., Fortran style)
      return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<2>,
                         nb::f_contig>(ptr, {rows, cols}, owner);
    },
    nb::arg("A"), nb::rv_policy::reference);
}

} // namespace sundials4py