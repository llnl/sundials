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

.. _Testing.GitLab:

GitLab CI Testing
=================

This section provides an overview of the GitLab continuous integration (CI)
features and describes the SUNDIALS GitLab CI setup. For more information see
the official `GitLab CI documentation <https://docs.gitlab.com/ee/ci/>`_.

For information specific to the LLNL GitLab CI see:

* `LLNL GitLab CI <https://lc.llnl.gov/confluence/display/GITLAB/GitLab+CI>`_

* `LLNL GitLab Runner Tags <https://lc.llnl.gov/gitlab/public-info/gitlab-ci/-/wikis/Gitlab-CI-Basic-Information>`_


SUNDIALS utilizes the GitLab CI pipeline code repository shared by LLNL RADIUSS
projects. The `docs <https://radiuss-shared-ci.readthedocs.io/en/latest/>`__ for
the shared project should be reviewed before reading the docs below.


CI Pipelines and Jobs
---------------------

Pipelines are the main component of the GitLab CI system and are comprised of
two components, stages and jobs. Stages define the logical steps in a pipeline
i.e., *when* to do something. Jobs are the basic configuration component and
define *what* to do. Each job is associated with a stage and multiple jobs in
the same stage are executed in parallel (if there are enough resources). If all
the jobs in a stage succeed, the pipeline moves on to the next stage. If any job
in a stage fails, the next stage is usually (see below) not executed and the
pipeline ends early.

Some pipelines are run automatically on new commits (after they are mirrored
from GitHub to LC GitLab). Other pipelines, such as the benchmarking pipeline,
are run on a schedule that is configured through the
`GitLab UI <https://lc.llnl.gov/gitlab/sundials/sundials/-/pipeline_schedules>`__.

Structure
^^^^^^^^^

As previously stated, most of our code for the LC GitLab CI pipelines is sourced from
the templates provided in the
`radiuss-shared-ci <https://radiuss-shared-ci.readthedocs.io/en/latest/>`__ repo.
Here we briefly outline the relevant files:

The ``.gitlab-ci.yml`` file in the root of the repository is the starting point for
defining the SUNDIALS GitLab CI pipelines. The only thing that is typically changed
in this file is the ``SPACK_REF`` variable in the ``variables`` section (this
is done when we update the version of Spack we use for installing dependencies).
Currently, we also override the `.build-and-test` job defined in this file so
that we can pull in some files from our fork of `radiuss-shared-ci`
(maintained `here <https://lc.llnl.gov/gitlab/sundials/radiuss-shared-ci>`__)
instead of the upstream repository.

The ``.gitlab/subscribed-pipelines.yml`` defines which machines we will test on.
This file may be modified if you need to add a new machine to test on.

The ``.gitlab/custom-jobs-and-variables.yml`` defines variables available in all
pipelines and jobs. This file may be modified if you need to add a new variable
that needs to be accessible to all pipelines and jobs.

The ``.gitlab/jobs/<machine>.yaml`` files define the jobs for a specific machine.
A "hidden" job of the form `.sundials_job_on_<machine>` is defined first in these
files and typically defines variables specific to that machine, such as what command
to use for executing MPI programs. The rest of the jobs in the file extend the
`.sundials_job_on_<machine>` and define the Spack spec that we will build and test.
Take for example, this Tioga job:

.. code-block:: YAML

   tioga_rocmcc_620_tpls:
     parallel:
       matrix:
         - COMPILER_SPEC: rocmcc@6.2.0
           AMDGPU_TARGET: [gfx90a]
     variables:
       # have to use ginkgo@master because our spack version does not have ginkgo@1.8.0: yet (which seems to be needed)
       # similarly, we need a newer magma than available to compile with 'rocm@6:' so we turn it off
       SPEC: "%${COMPILER_SPEC} cstd=99 cxxstd=17 precision=double amdgpu_target=${AMDGPU_TARGET} +rocm+mpi~magma+ginkgo+kokkos ^ginkgo@master+rocm amdgpu_target=${AMDGPU_TARGET} ^kokkos+rocm amdgpu_target=${AMDGPU_TARGET}"
     before_script:
       - module load rocmcc/6.2.0-magic
     extends: [.sundials_job_on_tioga]

The ``parallel:`` and ``matrix:`` keywords could be used to enable creating multiple jobs
with different variable values for each instance of the job, e.g., one job using
``rocmcc@6.2.0`` and another using ``rocmcc@6.2.1``. However, right now they only create
a single job (hence why ``COMPILER_SPEC`` and ``AMDGPU_TARGET`` only have one value). These
variables values are then used to create an environment variable ``SPEC`` which is the Spack spec
used by ``build_and_test.sh`` (discussed below) to configure and build SUNDIALS and the
necessary dependencies.

Disabling a Job
^^^^^^^^^^^^^^^

A job can be disabled by adding the variable ``.ON_<machine>: "OFF"`` to the ``variables:``
section of the job, e.g.,


.. code-block:: YAML

   tioga_rocmcc_620_tpls:
     parallel:
       matrix:
         - COMPILER_SPEC: rocmcc@6.2.0
           AMDGPU_TARGET: [gfx90a]
     variables:
       ON_TIOGA: "OFF" # disable this job
       SPEC: "%${COMPILER_SPEC} cstd=99 cxxstd=17 precision=double amdgpu_target=${AMDGPU_TARGET} +rocm+mpi~magma+ginkgo+kokkos ^ginkgo@master+rocm amdgpu_target=${AMDGPU_TARGET} ^kokkos+rocm amdgpu_target=${AMDGPU_TARGET}"
     before_script:
       - module load rocmcc/6.2.0-magic
     extends: [.sundials_job_on_tioga]

These variables can also be set when manually or scheduling a pipeline in the GitLab UI.

Updating Spack
^^^^^^^^^^^^^^

To update the spack commit used for the CI simply replace the commit hash in the
``SPACK_REF`` variable inside the ``.gitlab-ci.yml`` file with the new commit hash.
The first time a pipeline runs with a new ``SPACK_REF`` the pipeline will take longer than
normal as a new Spack build cache must be created and populated (so all packages will be
built from source).


Benchmark Jobs
^^^^^^^^^^^^^^

See :ref:`SUNDIALS Continuous Performance Testing (CPT)<CPT>` for more details.


Directories and Permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``/usr/workspace/sundials`` is the workspace directory associated with the ``sundials`` LC group.
  Users must be added to this group through the LC IDM application.

* ``/usr/workspace/sundials/ci`` is where all GitLab CI related files are stored.
  The correct permissions for this directory are ``drwxrws---``.

* ``/usr/workspace/sundials/ci/.builds`` is where GitLab CI pipelines are run. The permissions
  for this directory are ``drwxrwx---``, but directories within it must be ``drwx------``.
  Files within it should be ``-rw-rw----`` (can add ``x`` for group and owner as appropriate).

* ``/usr/workspace/sundials/ci/spack_stuff`` contains the Spack build caches amongst other Spack
  files. The permissions for this directory and directories below should be ``drwxrws---``.


GitLab CI Test Script
---------------------

The GitLab CI uses the script ``.gitlab/build_and_test.sh``, and when
benchmarking ``.gitlab/build_and_bench.sh``, to configure,
build, and test SUNDIALS. This script leverages two Git submodules:

* `uberenv <https://github.com/LLNL/uberenv>`_ -- automates using a package
  manager (e.g., Spack) to configure and build software. The top-level file
  ``.uberenv_config.json`` defines information need by uberenv including the
  the Spack commit to utilize and the location of Spack config and package
  files.

* `radiuss-spack-configs <https://github.com/sundials-codes/radiuss-spack-configs.git>`_
  -- is the SUNDIALS fork of the `LLNL radiuss-spack-configs <https://github.com/LLNL/radiuss-spack-configs>`_
  repository that provides spack environment files for various LLNL platfornms
  i.e., ``spack.yaml`` for Dane, Tioga, etc.

These submodules work in conjunction with ``scripts/sundials/package.py``
to configure and build any third-party libraries needed by the SUNDIALS
configuration and generates an initial CMake cache file for building SUNDIALS.
Other packages can be added to ``spack/packages`` if the default Spack package
needs to be overridden.

Spack Build Cache
^^^^^^^^^^^^^^^^^

The ``build_and_test.sh`` script leverage Spack build caches in ``/usr/workspace/sundials/ci/spack_stuff/<SPACK_REF>``
to speedup builds. These caches store binaries of packages that have been built previously. Separate caches are
made for each ``SPACK_REF`` to avoid conflicts across Spack versions.


Running Locally
^^^^^^^^^^^^^^^

It is possible to run these scripts locally on an LC machine. First set a ``SPACK_REF``
environment variable to a spack commit that you want to use, and then set a ``SPEC``
environment variable with a SUNDIALS spack spec that you want to test.
