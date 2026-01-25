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

.. _Testing:

Testing
=======

We include several types of tests in SUNDIALS: unit tests, integration test and performance tests.
These tests are run via as part of our Continuous Integration (CI) suite.

**Unit Tests**

The unit tests reside in two places: ``test/unit_tests`` and in files named ``test_`` within
``examples/``. With few exceptions, SUNDIALS unit tests should return an exit code indicating if the
test was successful or not. I.e., they should be self-contained and not require any output file to
determine success or failure.

**Integration Tests**

The integration tests are dual purpose; they serve as tests and also examples for users. They are
found in ``examples/``. The integration tests produce output that is checked against an output file,
i.e. an answer file, to determine if the test passed or failed. See :ref:`Answers <Answers>` for details.

**Performance Tests**

These tests are benchmarks of SUNDIALS performance and are found in ``benchmarks/``. Refer to
:ref:`Continuous Performance Testing <CPT>` for more detail.


.. toctree::
   :maxdepth: 1

   CTest
   Scripts
   GitHub
   GitLab
   Answers
   Benchmarking
   Spot
