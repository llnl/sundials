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
#include <sunmatrix/sunmatrix_sparse.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunmatrix_sparse(nb::module_& m)
{
#include "sunmatrix_sparse_generated.hpp"

  m.def(
    "SUNSparseMatrix_Data",
    [](SUNMatrix A)
    {
      auto nnz   = static_cast<size_t>(SUNSparseMatrix_NNZ(A));
      auto owner = nb::find(A);
      auto ptr   = SUNSparseMatrix_Data(A);
      // SUNSparseMatrix_Data returns data that cannot be directly indexed as a 2-dimensional numpy array
      return nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig>(ptr,
                                                                            {nnz},
                                                                            owner);
    },
    nb::arg("A"), nb::rv_policy::reference);

  m.def(
    "SUNSparseMatrix_IndexValues",
    [](SUNMatrix A)
    {
      auto nnz   = static_cast<size_t>(SUNSparseMatrix_NNZ(A));
      auto owner = nb::find(A);
      auto ptr   = SUNSparseMatrix_IndexValues(A);
      return nb::ndarray<sunindextype, nb::numpy, nb::ndim<1>, nb::c_contig>(ptr,
                                                                             {nnz},
                                                                             owner);
    },
    nb::arg("A"), nb::rv_policy::reference);

  m.def(
    "SUNSparseMatrix_IndexPointers",
    [](SUNMatrix A)
    {
      auto np    = static_cast<size_t>(SUNSparseMatrix_NP(A));
      auto owner = nb::find(A);
      auto ptr   = SUNSparseMatrix_IndexPointers(A);
      return nb::ndarray<sunindextype, nb::numpy, nb::ndim<1>,
                         nb::c_contig>(ptr, {np + 1}, owner);
    },
    nb::arg("A"), nb::rv_policy::reference);
}

} // namespace sundials4py
