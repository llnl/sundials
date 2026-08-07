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

.. _Test.GitHub:

GitHub CI Testing
=================

There are two categories of CI testing that we run on GitHub via `GitHub actions <https://github.com/LLNL/sundials/actions>`_:

1. Comprehensive (excluding GPUs)
2. Minimal/Short

The comprehensive testing is run only on one platform (Ubuntu) and utilizes Docker + Spack + the
latest E4S release to build SUNDIALS in many different configurations with and without third-party
libraries enabled.

The minimal/short testing runs on more platforms: Windows (MinGW and MSVC), and MacOS but it runs in
only one configuration (using 64-bit indices and double precision) and without third-party
libraries.


Building the Docker containers for CI
-------------------------------------

Original Setup
^^^^^^^^^^^^^^

These are the steps that were originally performed by Cody Balos
to build the Docker container(s) used for the comprehensive CI testing.

1. Create a ``spack.yaml`` for each configuration of SUNDIALS. E.g., int64 and double:

.. code-block:: yaml

  # This is a Spack Environment file.
  #
  # It describes a set of packages to be installed, along with
  # configuration settings.
  spack:
    packages:
      all:
        providers:
          blas: [openblas]
          mpi: [openmpi]
    # add package specs to the ``specs`` list
    specs:
    - hypre+int64~internal-superlu
    - petsc+double+int64
    - openmpi
    - openblas+ilp64
    - suite-sparse
    - superlu-dist+int64 ^parmetis+int64
    - trilinos+tpetra gotype=long_long
    config: {}
    modules:
      enable: []
    repos: []
    upstreams: {}
    container:
      images:
        os: ubuntu:20.04

2. Run ``spack containerize > Dockerfile`` in the directory of the ``spack.yaml``

3. The Dockerfile produced in step 2 was then manually modified to leverage
   Spack's buildcache feature and Docker's 'cache' bind-mount. The gist is that
   if the Spack build fails while building the Docker image, the buildcache
   makes it possible to reuse the binaries for the packages that were already installed
   before the build failed. Without the buildcache, the spack build failing would
   result in all packages needing to be built again when re-attempting to build the Docker image.

4. Run ``DOCKER_BUILDKIT docker build -t sundials-ci/<index-size>-<precision>:<tag>``

5. Push

Automated building of new containers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We currently have six different containers, one for each combination of {int32, int64} and {single,
double, extended} precision. These containers are pinned to an E4S release. When E4S does a release,
we can rebuild these containers to use the packages from it. We add E4S as a mirror in the Spack
environment so that its buildcache can be leveraged.

We also maintain two containers for the {int32, double} pair that are built automatically (in a
GitHub action) every week against the latest Spack develop commit. This allows us to test against
the latest versions of dependencies regularly and detect interface breakages.

Running Locally
---------------

It is possible to use the SUNDIALS CI containers available on `GitHub
<https://github.com/orgs/LLNL/packages?ecosystem=container>`__ for testing
locally. This will allow you to get as close as possible to running tests in the
CI environment on your own machine.

If you have Docker or `Podman <https://podman.io/>`__ installed on your machine already,
than the easiest way to use the containers is via the CMake targets ``setup_local_ci``
and ``test_local_ci``, for example:

.. code-block:: shell

   $ cd builddir
   $ make setup_local_ci
   ...
   $ make test_local_ci
   ...

The ``setup_local_ci`` target will pull the container for the configured ``SUNDIALS_PRECISION``
and ``SUNDIALS_INDEX_SIZE``. Then the ``test_local_ci`` target will run the full test suite
in the container.

The following CMake options are available for setting up local testing with the
containers.

.. cmakeoption:: SUNDIALS_TEST_CONTAINER_EXE

   Path to docker or podman

   Default: CMake will attempt to locate docker or podman automatically

.. cmakeoption:: SUNDIALS_TEST_CONTAINER_RUN_EXTRA_ARGS

   Extra arguments to pass to docker/podman run command

   Default: ``--tls-verify=false``

.. cmakeoption:: SUNDIALS_TEST_CONTAINER_MNT

   Path to project root inside the container

   Default: ``/sundials``

Alternatively, if you want to work with Docker directly, you can pull the image(s) and then
run the test suite manually. The ``run`` command will pull the image and start the container:

.. code-block:: shell

   docker run -t -d --name sundialsci-int32-double-latest --cap-add SYS_PTRACE -v "/path/to/your/sundials/development/repo":/sundials ghcr.io/llnl/sundials-ci-int32-double:latest

The ``exec`` command can then be used to execute the test script:

.. code-block:: shell

   docker exec -w /sundials/test sundialsci-int32-double-latest ./test_driver.sh --testtype CUSTOM --env env/docker.sh --tpls --sunrealtype double --indexsize 32

Alternatively, you can drop into a bash shell inside the container to run specific examples:

.. code-block:: shell

   docker exec -it sundialsci-int32-double-latest bash

On Macs, it is recommended to use `Podman <https://podman.io/>`__ (and then the
same steps above apply using ``podman`` instead of ``docker``). Podman is
useful on Linux too, as it can run rootless easily.

.. note::

   Its possible that the difference in your local machine architecture, and the
   one used to build the docker container(s), results in different answers and
   failing tests. You can provide the path to your own directory with answer
   files by setting the environment variable ``SUNDIALS_TEST_ANSWER_DIR`` with
   the path, and adding the argument ``-e SUNDIALS_TEST_ANSWER_DIR`` to the
   ``docker run`` command above.
