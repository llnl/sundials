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

.. _Test.Scripts:

CI Testing Scripts
==================

Some of the CI suites utilize various scripts in the ``test`` directory to setup
and run tests.

Test Driver
-----------

The script ``test_driver.sh`` is used to drive tests in the Jenkins and
GitHub CI. The input option ``--testtype`` can be used to select from a
predefined set of test configurations:

* BRANCH -- C99 build tests and a small set of build and run tests.

* PR -- C99 build tests and a larger set of build, run, and install tests using
  the SUNDIALS release tarball.

* RELEASE -- C99 build tests and an even larger set of build, run, and install
  tests using the SUNDIALS release tarball.

* CUSTOM -- run a single user-defined test configuration set using additional
  input options.

Run ``test_driver.sh -h`` for more information on the options available.

Testing Environment
-------------------

The environment setup script, ``env/setup_env.sh``, and a machine-specific
script are used by the test driver to configure the testing environment. The
machine script can be specified using the ``--env`` option to ``test_driver.sh``
otherwise, the following variables/locations will be checked:

#. The script specified by the environment variable ``SUNDIALS_ENV_FILE``
#. A user's local environment script: ``<sunrepo>/test/env/env.sh``

The setup script will set various environment variables based on inputs from the
test driver to configure SUNDIALS. Any unrecognized input options passed to the
test driver will be passed through to the machine script to supply additional
setup information. For example, a compiler spec could be passed through to the
machine script to change the compiler used when testing.

The setup script will set the following environment variables (among others):

* ``SUNDIALS_PRECISION`` the real type: ``single``, ``double``, or ``extended``

* ``SUNDIALS_INDEX_SIZE`` the index size: ``32`` or ``64``

* ``SUNDIALS_LIBRARY_TYPE`` the library type: ``shared``, ``static`` or ``both``

* ``SUNDIALS_TPLS`` enable third-party libraries (TPLs): ``ON`` or ``OFF``

* ``SUNDIALS_TEST_TYPE`` enable additional tests: ``STD`` or ``DEV``

  - When set to ``DEV``, development and unit tests will be enabled

After this initial setup, the machine script is called to set any additional
environment variables need to configure SUNDIALS for testing. Currently, these
include the following:

.. code-block:: shell

   CC  = C compiler
   CXX = C++ compiler
   FC  = Fortran compiler

   CFLAGS   = C compiler flags
   CXXFLAGS = C++ compiler flags
   FFLAGS   = Fortran compiler flags

   # Enable/disable MPI support
   SUNDIALS_MPI = ON or OFF
   MPICC        = MPI C compiler wrapper
   MPICXX       = MPI C++ compiler wrapper
   MPIFC        = MPI Fortran compiler wrapper
   MPIEXEC      = executable for launching MPI runs

   # Enable/disable PThread support
   SUNDIALS_PTHREAD = ON or OFF

   # Enable/disable OpenMP support
   SUNDIALS_OPENMP = ON or OFF

   # Enable/disable OpenMP device offloading support
   SUNDIALS_OPENMPDEV = ON or OFF

   # Enable/disable CUDA support
   SUNDIALS_CUDA = ON or OFF

   # Enable/disable HIP support
   SUNDIALS_HIP = ON or OFF

   # Enable/disable RAJA support
   SUNDIALS_RAJA = ON or OFF
   RAJA_ROOT     = full path to RAJA installation
   RAJA_BACKENDS = RAJA backends

   # Enable/disable SYCL support
   SUNDIALS_SYCL = ON or OFF

   # Enable/disable LAPACK linear solvers
   SUNDIALS_LAPACK  = ON or OFF
   LAPACK_LIBRARIES = full path to LAPACK library

   # Enable/disable KLU linear solver
   SUNDIALS_KLU             = ON or OFF
   SUITE_SPARSE_INCLUDE_DIR = full path to SuiteSparse include directory
   SUITE_SPARSE_LIBRARY_DIR = full path to SuiteSparse library directory

   # Enable/disable SuperLU_MT linear solver
   SUNDIALS_SUPERLU_MT    = ON or OFF
   SUPERLU_MT_INCLUDE_DIR = full path to SuperLU_MT include directory
   SUPERLU_MT_LIBRARY_DIR = full path to SuperLU_MT library directory

   # Enable/disable SuperLU_DIST linear solver
   SUNDIALS_SUPERLU_DIST  = ON or OFF
   SUPERLU_DIST_INCLUDE_DIR = full path to SuperLU_DIST include directory
   SUPERLU_DIST_LIBRARY_DIR = full path to SuperLU_DIST library directory
   SUPERLU_DIST_LIBRARIES   = additional link libraries for SuperLU_DIST

   # Enable/disable MAGMA linear solver
   SUNDIALS_MAGMA = ON or OFF
   MAGMA_ROOT     = full path to MAGMA installation
   MAGMA_BAKCENDS = MAGMA backend

   # Enable/disable hypre support
   SUNDIALS_HYPRE    = ON or OFF
   HYPRE_INCLUDE_DIR = full path to hypre include directory
   HYPRE_LIBRARY_DIR = full path to hypre library directory

   # Enable/disable PETSc support
   SUNDIALS_PETSC = ON or OFF
   PETSC_ROOT     = full path to PETSc installation

   # Enable/disable Trilinos support
   SUNDIALS_TRILINOS = ON or OFF
   TRILINOS_ROOT     = full path to Trilinos installation

   # Enable/disable Trilinos support
   SUNDIALS_XBRAID = ON or OFF
   XBRAID_ROOT     = full path to XBraid installation

Test Runner
-----------

When comparing test outputs with answer files or profiling tests with Caliper
the Python test runner, ``test/testRunner``, is used run the test under CTest.
