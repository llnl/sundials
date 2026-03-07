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

.. _Benchmarks:

Benchmarks
==========

In this chapter we describe performance benchmark programs distributed with
SUNDIALS. To enable the benchmark programs and set the installation location use
the following CMake options.

.. cmakeoption:: SUNDIALS_ENABLE_BENCHMARKS

   Enable the benchmark problems

   Default: ``OFF``

   .. versionadded:: x.y.z

      Replaces the deprecated option ``BUILD_BENCHMARKS``

.. cmakeoption:: SUNDIALS_BENCHMARKS_INSTALL_PATH

   Full path to where to install the benchmark problems

   Default: ``CMAKE_INSTALL_PREFIX/benchmarks``

   .. versionadded:: x.y.z

      Replaces the deprecated option ``BENCHMARKS_INSTALL_PATH``

The following sections describe the SUNDIALS performance benchmark problems.

.. toctree::
   :maxdepth: 1

   advection_reaction.rst
   diffusion.rst
