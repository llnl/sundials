..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. warning::

   sundials4py is in beta and is subject to breaking changes. We welcome
   feedback on this new feature of SUNDIALS. Please report issues via the
   :github:`SUNDIALS GitHub repository <issues>`.

`sundials4py <https://pypi.org/project/sundials4py/>`__ provides official
(supported by the SUNDIALS team) Python bindings to much of the SUNDIALS
library, allowing you to use SUNDIALS directly from Python.

The bindings are built using `nanobind <https://nanobind.readthedocs.io/>`__ and
`litgen <https://pthom.github.io/litgen>`__ and are designed to be easy to use
from Python in conjunction with ubiquitous libraries in the Python scientific
computing and machine learning ecosystems. To that end, sundials4py supports:

- Python's automatic memory management
- Python definitions of user-supplied callback functions
- Zero-copy exchange of arrays (CPU and Device) through DLPack protocol and
  numpy's ndarray

.. note::

   sundials4py requires Python 3.12+

The Python User Guide focuses on specific aspects of using SUNDIALS from Python
and assumes the user is familiar with SUNDIALS. New SUNDIALS users should first
read the :ref:`General User Guide <SUNDIALS>` to understand the features and
usage of SUNDIALS packages.
