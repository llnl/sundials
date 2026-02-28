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

SPOT Performance Visualization Tool
====================================

In order to protect against performance regression of SUNDIALS at all scales,
we leverage the GitLab CI infrastructure setup for CI testing to perform
continuous performance testing.

To track and visualize performance over time, we collect data with
`Caliper <https://lc.llnl.gov/confluence/display/CALI/Spot+DB>`_ and
`Adiak <https://github.com/LLNL/Adiak>`_. Steps to setting up continuous
performance testing(CPT) can be found in
:ref:`Continuous Performance Testing<CPT>`.

Documentation to using the SPOT performance visualization tool can be found
`here <https://lc.llnl.gov/confluence/display/SpotDoc/Spot+Documentation>`_.

Setting Up Your Own SPOT Visualizations
---------------------------------------

To display data from non-GitLab jobs or create a local collection of runs, input
into the SPOT search bar the absolute path to the directory containing the
``.cali`` files and refresh the page. By default ``.cali`` files will be
generated in the build directory under ``Benchmarking/caliper`` and
``Testing/caliper``. To specify where ``.cali`` output files are placed, define
the CMake options ``SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR`` and
``SUNDIALS_TEST_CALIPER_OUTPUT_DIR`` with the desired directory path.

To retain the same filters as a given SPOT visualization link in
:ref:`Bookmarks`, swap out the ``sf`` value in the URL with the
directory path containing the ``.cali`` files.

.. _Bookmarks:

Bookmarks to Notable SPOT Visualizations
----------------------------------------------

Notes:

SPOT will be put into maintenance mode in the coming months
(as of July 2023). The team behind SPOT plans to bring the functionality to
`Thicket <https://github.com/llnl/thicket>`_ as SPOT's successor.

Aggregate results for example runs on SPOT may not be accurate. The aggregation method
involves combining ``.cali`` files in ways that do not necessarily equal the
true sum of the results.

**Examples**

`Aggregated Runs <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Example&ch_executable=1&ch_launchdate=1&groupby=cluster&aggregate=avg&xaxis=job_start_time&yaxis=Max%20time%2Frank>`_

**Benchmarks**

**Advection Reaction 3D**

* `All configurations <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D&ch_executable=1&ch_launchdate=1&groupby=cmdline&yaxis=Max%20time%2Frank&aggregate=avg>`__

* By configuration

  * Serial

    * `ARK-DIRK, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_arkdirk_newton&ch_executable=1&ch_launchdate=1&aggregate=&yaxis=Max%20time%2Frank&groupby=cluster>`__

    * `ARK-IMEX, TL-Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_arkimex_tlnewton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank&groupby=cluster>`__

    * `CV-BDF, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_cvbdf_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank&groupby=cluster>`__

    * `IDA, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_ida_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank&groupby=cluster>`__

  * CUDA

    * `ARK-DIRK, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpicuda_arkdirk_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

    * `ARK-IMEX, TL-Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpicuda_arkimex_tlnewton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

    * `CV-BDF, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpicuda_cvbdf_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

    * `IDA, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpicuda_ida_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

  * HIP

    * `ARK-DIRK, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpihip_arkdirk_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

    * `ARK-IMEX, TL-Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpihip_arkimex_tlnewton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

    * `CV-BDF, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpihip_cvbdf_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

    * `IDA, Newton <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/advection_reaction_3D/advection_reaction_3D_raja_mpihip_ida_newton&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__


**Diffusion 2D**

Note: CUDA Diffusion 2D visualizations are not available as the benchmark errors out before completion.

* `All configurations <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D&ch_executable=1&ch_launchdate=1&groupby=executable&yaxis=Max%20time%2Frank&aggregate=avg>`__

* MPI + Serial

  * `ARKODE <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D/arkode_diffusion_2D_mpi_d2d_arkode_serial&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

  * `CVODE <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D/cvode_diffusion_2D_mpi_d2d_cvode_serial&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

  * `IDA <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D/ida_diffusion_2D_mpi_d2d_ida_serial&ch_executable=1&ch_launchdate=1&yaxis=Max%20time%2Frank>`__

* MPI + HIP

  * `ARKODE <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D/arkode_diffusion_2D_mpihip_d2d_arkode_hip&ch_executable=1&ch_launchdate=1&aggregate=max>`__

  * `CVODE <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D/cvode_diffusion_2D_mpihip_d2d_cvode_hip&ch_executable=1&ch_launchdate=1&aggregate=max>`__

  * `IDA <https://lc.llnl.gov/spot2/?sf=/usr/workspace/sundials/califiles/Benchmarking/diffusion_2D/ida_diffusion_2D_mpihip_d2d_ida_hip&ch_executable=1&ch_launchdate=1&aggregate=max>`__
