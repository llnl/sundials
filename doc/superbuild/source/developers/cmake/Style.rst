..
   Author(s): David J. Gardner @ LLNL
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

.. _CMake.Style:

Style
=====

In this section we describe the CMake style conventions and guidelines for
SUNDIALS.

Naming
------

SUNDIALS CMake cache variables (configuration options) must adhere to the
following naming conventions:

#. To avoid naming collisions with other projects that include SUNDIALS as part
   of their CMake build system, SUNDIALS CMake cache variables are prefixed with
   ``SUNDIALS_``.

#. Cache variables use screaming snake case, ``SUNDIALS_OPTION_NAME``.

#. The ``SUNDIALS_`` prefix may be followed by a category name for grouping
   related options, ``SUNDIALS_CATEGORY_OPTION_NAME``. For example,

   * ``SUNDIALS_TEST_`` for SUNDIALS testing options

   * ``SUNDIALS_EXAMPLES_`` for options related to the SUNDIALS examples

   * ``SUNDIALS_BENCHMARKS_`` for options related to SUNDIALS benchmarks

   * ``SUNDIALS_DEV_`` for development options

   * ``SUNDIALS_DEBUG_`` for debugging options

#. Boolean options are named ``SUNDIALS_ENABLE_OPTION_NAME`` or
   ``SUNDIALS_CATEGORY_ENABLE_OPTION_NAME``.

#. Options to enable third party libraries are named
   ``SUNDIALS_ENABLE_LIBRARY_NAME``.

Formatting
----------

CMake files should be formatted with `cmake-format
<https://cmake-format.readthedocs.io>`__. ``cmake-format`` can be installed with
``pip``:

.. code-block:: console

   pip install cmake-format

Individual files can be formatted with the command:

.. code-block:: console

   cmake-format -i <cmake file>

The ``.cmake-format.py`` file in the root of the project defines our
configuration for ``cmake-format``.
