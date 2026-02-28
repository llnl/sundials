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

.. _OpenPR:

Opening a Pull Request
======================

When a branch is ready to be reviewed for integration into the ``main`` or
``develop`` branches follow the steps below to open a pull request:

#. Browse to `https://github.com/LLNL/sundials <https://github.com/LLNL/sundials>`_

#. Click on the branch icon on the left side of screen - you'll see a list
   of all available branches

#. Click on your branch - you'll see a 'Compare' screen that lets you pick a
   branch (source on top) to merge with another branch (target on bottom)

#. Select the desired branches and click 'Create pull request'

#. Edit the title of the pull request (defaults to branch name), add a
   description, and select reviewers that can approve the request

#. Click 'Create'

The selected reviewers will go over the changes in the pull request and may
ask for additional changes before merging the branch. After the pull request is
merged, delete the local copy the branch:

.. code-block:: none

   $ git checkout PARENT
   $ git branch -D <branch-name>


.. note::

   When you open a PR, various "checks" will run on your code. You can see the status of these
   checks at the bottom of the PR page on GitHub. Any fixes needed to pass the checks for
   formatting, spelling, and SWIG updates can be downloaded as git patches under "artifacts".
   These patches can be applied to your branch using ``git apply <patch>`` and then committed.

.. warning::

   Almost all pull requests should be issued against the ``develop`` branch.
   Generally, the only branches we merge to ``main`` are special PRs to synchronize
   ``develop`` and ``main``.
