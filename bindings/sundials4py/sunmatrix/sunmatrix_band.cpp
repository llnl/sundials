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
#include <sunmatrix/sunmatrix_band.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunmatrix_band(nb::module_& m)
{
#include "sunmatrix_band_generated.hpp"

  m.def(
    "SUNBandMatrix_Data",
    [](SUNMatrix A)
    {
      auto ldata = static_cast<size_t>(SUNBandMatrix_LData(A));
      auto owner = nb::find(A);
      auto ptr   = SUNBandMatrix_Data(A);
      // SUNBandMatrix_Data returns data that cannot be directly indexed as a 2-dimensional numpy array
      return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig>(ptr,
                                                                            {ldata},
                                                                            owner);
    },
    nb::arg("A"), nb::rv_policy::reference);
}

} // namespace sundials4py
