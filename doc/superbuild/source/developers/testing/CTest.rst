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

.. _Testing.CTest:

CTest
=====

We use `CTest <https://cmake.org/cmake/help/latest/manual/ctest.1.html>`__ to
run tests with the CMake ``test`` target. For example, the default testing
configuration can be run with:

.. code-block:: shell

   mkdir build
   cd build
   cmake ../.
   make
   make test

In this case a subset of the example programs (if examples are enabled) are run
as a smoke test. Most examples pass if they run successfully while examples with
an analytic solution may check the solution error and return pass or fail
accordingly.

Development and Unit Tests
--------------------------

Beyond the basic smoke tests, development and unit tests can be enabled with
following CMake options.

.. cmakeoption:: SUNDIALS_TEST_ENABLE_DEV_TESTS

   Include development tests

   Default: OFF

.. cmakeoption:: SUNDIALS_TEST_ENABLE_UNIT_TESTS

   Include unit tests

   Default: OFF

.. cmakeoption:: SUNDIALS_TEST_ENABLE_GTEST

   Include unit tests that utilize `Google test <https://github.com/google/googletest>`__

   When enabled, CMake will automatically download Google test.

   Default: ON (when unit tests are enabled)

Answer Files
------------

Generally, unit tests check for correctness internally and return pass or fail
accordingly. Most development tests, like the example programs, pass if they run
successfully. To check for regressions, we compare the output from tests and
examples that do not internally verify correctness with saved "answer" files
that are deemed correct by the SUNDIALS team. The following CMake options can be
used to enable and configure output comparisons.

.. cmakeoption:: SUNDIALS_TEST_ENABLE_DIFF_OUTPUT

   Compare test outputs with answer files

   Default: ON (when development or unit tests are enabled)

.. cmakeoption:: SUNDIALS_TEST_FLOAT_PRECISION

   Precision for floating point comparisons (number of digits)

   Default: 4

.. cmakeoption:: SUNDIALS_TEST_INTEGER_PRECISION

   Precision for integer comparisons (percent difference)

   Default: 10

Due to differences in hardware and floating point round-off, it is possible that
test outputs on one machine will differ from the answers generated on another
leading to erroneous failures. As such, the following CMake option can be used
to specify a directory containing alternative answer files generated on the same
machine from the ``develop`` branch.

.. cmakeoption:: SUNDIALS_TEST_ANSWER_DIR

   Location of test answer files

   Default: Use output (``.out``) files in the same directory as the test source

To assist in creating answer files for a new machine, the CMake option below can
be used to change the directory where test output files are written when running
tests.

.. cmakeoption:: SUNDIALS_TEST_OUTPUT_DIR

   Location to write test output files

   Default: ``<cmake build directory>/Testing/output``

For example, answer files for a minimal configuration can be generated with the
following steps.

.. code-block:: shell

   git checkout develop
   mkdir build
   cd build
   cmake ../. \
     -DSUNDIALS_TEST_ENABLE_DEV_TESTS=ON \
     -DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON \
     -DSUNDIALS_TEST_OUTPUT_DIR=<machine output directory>
   make
   make test

Compiler Flags
--------------

The CI suites discussed later build SUNDIALS with additional compiler warnings
enabled using the following CMake options. The specific flags used depend on the
real type precision, index size, and if the Fortran interfaces are enabled. See
``cmake/SundialsSetupCompilers.cmake`` for the exact set of flags. The current
set of flags is compatible with with GNU and Clang compilers.

.. cmakeoption:: ENABLE_ALL_WARNINGS

   Enable additional compiler warnings

   Default: OFF

.. cmakeoption:: ENABLE_WARNINGS_AS_ERRORS

   Treat compiler warnings as errors

   Default: OFF

Additionally, the CI will run a subset of tests (no TPLs) using different
sanitizers. These can be enabled with the following CMake options.

.. cmakeoption:: ENABLE_ADDRESS_SANITIZER

   Enable sanitizer to detect memory errors, adds the ``-fsanitize=address``
   flag. Depending on the compiler, this may also detect memory leaks.

   Default: OFF

.. cmakeoption:: ENABLE_LEAK_SANITIZER

   Enable sanitizer to detect memory leaks, adds the ``-fsanitize=leak``
   flag. Depending on the compiler, the leak sanitizer may be part of the
   address sanitizer.

   Default: OFF

.. cmakeoption:: ENABLE_MEMORY_SANITIZER

   Enable sanitizer to detect uninitialized memory errors, adds the
   ``-fsanitize=memory`` flag.

   Default: OFF

.. cmakeoption:: ENABLE_UNDEFINED_BEHAVIOR_SANITIZER

   Enable sanitizer to detect undefined behavior errors, adds the
   ``-fsanitize=undefined`` flag.

   Default: OFF
