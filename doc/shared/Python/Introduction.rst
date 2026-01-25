.. ----------------------------------------------------------------
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
   ----------------------------------------------------------------

.. _Python.Introduction:

Introduction
============

.. warning::

   sundials4py is currently in beta. The module API is subject to change even in SUNDIALS patch releases.


sundials4py is designed to be easy to use from Python in conjunction with ubiquitous libraries in the Python scientific computing and machine learning ecosystems.
To that end, it supports:

- Python's automatic memory management
- Python definitions of user-supplied callback functions
- Zero-copy exchange of arrays (CPU and Device) through DLPack protocol and numpy's ndarray

sundials4py is built using `nanobind <https://nanobind.readthedocs.io/en/latest/index.html#>`__ and `litgen <https://pthom.github.io/litgen/litgen_book/00_00_intro.html>`__.

.. note::

   sundials4py requires Python 3.12+


Installation
------------

You can install sundials4py directly from `PyPI <https://pypi.org/project/sundials4py/>`__ using pip:

.. code-block:: bash

   pip install sundials4py

Or, you can install from git:

.. code-block:: bash

   pip install git+https://github.com/LLNL/sundials.git

The default build of sundials4py that is distributed as a binary wheel uses double precision real types and 64-bit indices. 
To install SUNDIALS with different precisions and index sizes, you can build from source wheels instead of using the pre-built 
binary wheels. When building from source wheels instead of binary wheels, you can customize the SUNDIALS precision (real type)
and index type at build time by passing the CMake arguments in an environment variable when running pip. For example:

.. code-block:: bash

   export CMAKE_ARGS="-DSUNDIALS_PRECISION=SINGLE -DSUNDIALS_INDEX_SIZE=64"
   pip install sundials4py --no-binary=sundials4py

Other SUNDIALS options can also be accessed in this way. Review :numref:`Installation.Options` for more information on the available options.

.. note::

   Not all SUNDIALS options are supported by the Python interfaces. In particular, third-party libraries are not yet supported.

After installation, you can import sundials4py with:

.. code-block:: python

   import sundials4py

sundials4py is comprised of the following submodules:

- ``sundials4py.core``: contains all the shared SUNDIALS classes and functions
- ``sundials4py.arkode``: contains all of the ARKODE specific classes and functions
- ``sundials4py.cvodes``: contains all of the CVODES specific classes and functions
- ``sundials4py.idas``: contains all of the IDAS specific classes and functions
- ``sundials4py.kinsol``: contains all of the KINSOL specific classes and functions

CVODE and IDA dot not have modules because CVODES and IDAS provide all of the same capabilities plus continuous forward and adjoint sensitivity analysis.

For more information on usage, differences from the C/C++ API and examples, continue to the next sections of this documentation.

.. _Python.Table:

.. list-table:: List of SUNDIALS Modules Available from Python
   :align: center
   :header-rows: 1
   :widths: 30 40

   * - **Class/Module**
     - **Python Module Name**
   * - SUNDIALS core
     - ``sundials4py.core``
   * - ARKODE
     - ``sundials4py.arkode``
   * - ARKODE::ARKSTEP
     - ``sundials4py.arkode``
   * - ARKODE::ERKSTEP
     - ``sundials4py.arkode``
   * - ARKODE::MRISTEP
     - ``sundials4py.arkode``
   * - ARKODE::SPRKSTEP
     - ``sundials4py.arkode``
   * - ARKODE::LSRKSTEP
     - ``sundials4py.arkode``
   * - ARKODE::SPLITTINGSTEP
     - ``sundials4py.arkode``
   * - ARKODE::FORCINGSTEP
     - ``sundials4py.arkode``
   * - CVODE
     - ``sundials4py.cvode``
   * - CVODES
     - ``sundials4py.cvodes``
   * - IDA
     - ``sundials4py.ida``
   * - IDAS
     - ``sundials4py.idas``
   * - KINSOL
     - ``sundials4py.kinsol``
   * - NVECTOR_MANVECTOR
     - ``sundials4py.core``
   * - NVECTOR_SERIAL
     - ``sundials4py.core``
   * - SUNADAPTCONTROLLER_IMEXGUS
     - ``sundials4py.core``
   * - SUNADAPTCONTROLLER_SODERLIND
     - ``sundials4py.core``
   * - SUNADAPTCONTROLLER_MRIHTOL
     - ``sundials4py.core``
   * - SUNADJOINTCHECKPOINTSCHEME_FIXED
     - ``sundials4py.core``
   * - SUNDOMEIGEST_POWER
     - ``sundials4py.core``
   * - SUNLINSOL_BAND
     - ``sundials4py.core``
   * - SUNLINSOL_DENSE
     - ``sundials4py.core``
   * - SUNLINSOL_PCG
     - ``sundials4py.core``
   * - SUNLINSOL_SPBCGS
     - ``sundials4py.core``
   * - SUNLINSOL_SPFGMR
     - ``sundials4py.core``
   * - SUNLINSOL_SPGMR
     - ``sundials4py.core``
   * - SUNLINSOL_SPTFQMR
     - ``sundials4py.core``
   * - SUNMATRIX_BAND
     - ``sundials4py.core``
   * - SUNMATRIX_DENSE
     - ``sundials4py.core``
   * - SUNMATRIX_SPARSE
     - ``sundials4py.core``
   * - SUNNONLINSOL_FIXEDPOINT
     - ``sundials4py.core``
   * - SUNNONLINSOL_NEWTON
     - ``sundials4py.core``
