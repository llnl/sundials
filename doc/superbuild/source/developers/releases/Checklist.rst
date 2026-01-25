..
   Author(s): David J. Gardner, Cody J. Balos @ LLNL
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

.. _Versioning:

Versioning
==========

SUNDIALS follows the `semantic versioning <https://semver.org/>`__ scheme and
each release is given a version number ``x.y.z`` where

* ``x`` is the major version number and is incremented when the release includes
  an incompatible API change

* ``y`` is the minor version number and is incremented when new functionality
  has been added and is backwards compatible with the prior version

* ``z`` is the patch version and is incremented when the release only includes
  backwards compatible bug fixes

Additionally a pre-release versions will have ``-l.w`` appended to the version
number i.e., ``x.y.z-l.w`` where

* ``l`` is a pre-release label e.g., ``alpha``, ``beta``, or ``rc``.

* ``w`` is the pre-release version e.g, ``alpha.0``, ``alpha.1``, etc.

The convention for shared library versioning is ``libfoo.so.x.y.z``. Note that
there was no library versioning before sundials-2.2.0. The version information
is specified by setting the ``VERSION`` and ``SOVERSION`` properties to the
library target

.. code-block:: cmake

   add_library(foo file1.c file2.c)
   set_target_properties(foo PROPERTIES VERSION x.y.z SOVERSION x)

and build system is responsible for creating the symlinks

.. code-block:: none

   libfoo.so   -> libfoo.so.x.y.z
   libfoo.so.x -> libfoo.so.x.y.z

Note that we force ``SOVERSION`` to be the same as the major number for
compatibility with the libtool versioning scheme (at least under Linux).

.. _ReleaseCheckList:

Pre-Release Tasks
=================

This is a list of tasks that need to be done before a SUNDIALS release.
The order is bottom-up starting with solver source files and ending with
web pages.

#. Create the release branch ``release/vX.Y.Z`` from ``develop`` in the git
   repo.

   .. code-block:: shell

      git checkout develop
      git pull # develop should be up to date with origin/develop now
      git checkout -b release/vX.Y.Z

#. Update version numbers in ``scripts/updateVerson.sh``.

#. Run the ``updateVerson.sh`` script.

   .. code-block:: shell

      pushd scripts/ && ./updateVersion.sh && popd

#. If this is a major release, search the SUNDIALS code for 'DEPRECATION NOTICE'
   and 'SUNDIALS_DEPRECATED'. All deprecated functions should be removed (unless
   this is the first version that they are deprecated).

#. Remove any empty sections from the ``CHANGELOG.md`` and
   ``doc/shared/RecentChanges.rst`` files.

#. Update version numbers of third party libraries in the Install Guide in doc
   directory.

#. Open a pull request from the release branch to ``develop`` with the title
   "Release: vX.Y.Z" and description "SUNDIALS Release vX.Y.Z". Add the pull
   request to the "SUNDIALS Next" milestone.

Release Procedure
=================

#. Once the release PR is passing all tests and has been approved, merge
   it. Like all merges to ``develop``, this should be done with the "Squash and
   merge" option.

#. Sync the main branch with develop. This merge is done locally rather than
   through a GitHub PR (so that the merge is a fast-forward). The steps are as
   follows:

   .. code-block:: shell

      git checkout develop
      git pull # develop should be up to date with origin/develop now
      git checkout main
      git pull # main should be up to date with origin/main now
      git merge --ff-only develop # we want to do a fast-forward merge (no merge commit)
      git tag -a vX.Y.Z -m 'SUNDIALS Release vX.Y.Z'
      git push --tags origin main

   .. note::

      The final step (pushing to main) requires temporarily changing the GitHub
      repository settings to allow whoever is doing the push an exception in the
      ``main`` branch protection rules.

   .. danger::

      Remember to remove this exception to the branch protection rules after
      making the push to update ``main``.

#. Once readthedocs finishes building the `"latest" release documentation
   <https://app.readthedocs.org/projects/sundials/>`__, create the release
   tarballs *on a Linux machine*.

   .. code-block:: shell

      pushd scripts/ && ./tarscript.sh && popd

   .. warning::

      Creating the tarballs on a Mac can cause issues. Furthermore, it is
      important to wait to create the tarballs until readthedocs finishes
      building the latest release docs so that cross-references have valid
      links.

   The script will compile the documentation PDFs (user guides and example docs)
   and create all the tarballs in their final form, appropriate for uploading as
   artifacts to the GitHub release, in the top-level ``tarballs`` directory.

   .. note::

      If you get an error about downloading the Sphinx ``objects.inv`` file
      during the "superbuild" documentation build, you can manually download the
      file with the command

      .. code-block:: shell

         pushd doc/superbuild && wget https://www.sphinx-doc.org/en/master/objects.inv && popd

#. Draft the `release on GitHub <https://github.com/LLNL/sundials/releases>`__.
   Select the tag for ``vX.Y.Z`` and use the title "SUNDIALS vX.Y.Z". The
   description of the release is just a copy of the ``CHANGELOG.md`` notes for
   the release with hard line-wraps removed. Attach the tarballs as well as the
   example documentation PDFs as binary artifacts.

#. On the `GitHub milestones page
   <https://github.com/LLNL/sundials/milestones>`__ rename the "SUNDIALS Next"
   milestone "SUNDIALS X.Y.Z", close the renamed milestone, and create a new
   "SUNDIALS Next" milestone.

#. Open a pull request to the `spack-packages
   <https://github.com/spack/spack-packages>`__ repo to add the latest release
   to the list of versions in the `SUNDIALS package
   <https://github.com/spack/spack-packages/blob/develop/repos/spack_repo/builtin/packages/sundials/package.py>`__
   e.g.,

   .. code-block:: python

      version("X.Y.Z", tag="vX.Y.Z", commit="vX.Y.Z commit hash")

#. Update Internal Drupal Web pages for SUNDIALS:
   https://computing-staging.llnl.gov/user

   a) Modify content (save edits on each page as you go)

      * Edit Main Page:
        https://computing-staging.llnl.gov/projects/sundials

      * Edit Download Page:
        https://computing-staging.llnl.gov/projects/sundials/sundials-software

        * Update main download table with links to new versions of solvers.
        * The example documentation links need to be updated as well.
        * Update Previous releases table with new entry for previous release of full SUNDIALS suite.

   b) Once each sub page is complete, ask for team review of draft pages:
      https://computing-staging.llnl.gov/projects/sundials

   c) After team comments are included and saved, select the
      "Publishing options" button in the bottom left group of buttons on the
      draft page. Select the Moderation state reflecting the amount of
      required review then Save. This must be done for each page and is the
      final action before pages are uploaded for external release.

#. After final push, ensure web content and behavior is as expected on the main
   page: http://computing.llnl.gov/projects/sundials

Post-Release Tasks
==================

#. Prepare SUNDIALS for the next release cycle using the following steps:

   .. code-block:: shell

      git checkout develop
      git checkout -b maintenance/start-new-release-cycle
      pushd scripts/ && ./startReleaseCycle.sh && popd
      git add . && git commit -m 'start new release cycle'
      git push -u origin maintenance/start-new-release-cycle
      # Now open the PR to develop on GitHub.
