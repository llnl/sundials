..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Workflow:

Workflow
========

This section describes the typical SUNDIALS workflow for adding new features and
patching released code based on the `Gitflow workflow
<https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow>`_.
For details on setting up Git and cloning the repository see the :ref:`GitSetup`
section. A list of helpful commands for working with Git is available in the
:ref:`GitCheatSheet`.

.. _Gitflow:

Gitflow
-------

In the Gitflow workflow there are two persistent branches in the repository,
``master`` and ``develop``. The ``master`` branch contains the official release
history while ``develop`` is an integration branch for new features. Work is not
done directly on ``master`` or ``develop``, rather changes to either branch are
merged in through pull requests from short-lived branches that engage a *single*
feature or issue. These working branches generally fall into two groups:

* **Feature branches** for developing new features for future releases
  (typically created off of ``develop``).

* **Maintenance branches** for patching the current production release
  (typically created off of ``master``).

New feature or maintenance branches are created and worked on until they are
ready to be merged into ``master`` or ``develop``. To merge the code changes
into the desired branch a pull request is created and the code changes are
reviewed by other team members. Once the pull request is approved, the changes
will be merged and the working branch can be deleted. See
the :ref:`ExampleWorkflow` sections for an example of a typical workflow.

.. _NamingBranches:

Naming Branches
---------------

To keep track of what work is in progress on different branches it is helpful to
use short, descriptive branch names. When creating a new feature branch, use the
``feature/`` label followed by a name for the feature e.g.,
``feature/HypreVector`` for adding for a *hypre* vector. For maintenance
branches, the branch label is usually ``bugfix/`` or ``maintenance/`` followed
by the corresponding JIRA ticket name e.g., ``bugfix/SUNDIALS-123`` or a name
that describes the issue being addressed e.g.,
``maintenance/documentation-updates``.

.. _ExampleWorkflow:

Example Workflow
----------------

In the following example workflow PARENT refers the branch from which the
current working branch was created i.e., the ``develop`` branch for feature
branches and the ``master`` for maintenance branches.

#. To start working, checkout an existing branch or create a new one.

   a. To work on an existing branch:

      i. Checkout the desired branch:

          .. code-block:: none

             $ git checkout <branch-name>

      ii. Make sure the branch is up to date:

          .. code-block:: none

             $ git pull

   b. To create a new branch to work on:

      i. Checkout the desired PARENT branch and make sure it's up to date:

         .. code-block:: none

            $ git checkout PARENT
            $ git pull

      ii. Create a new branch off of the PARENT branch:

          .. code-block:: none

             $ git checkout -b <branch-name> PARENT

      iii. Push the new branch to the remote repo:

           .. code-block:: none

              $ git push -u origin <branch-name>

   NOTE: When checking out a branch, any local changed must be committed (or
   stashed) prior to switching branches.

#. Modify files with your editor of choice:

   .. code-block:: none

      $ emacs <some-path/file-name>

#. Stage the files you want to commit with:

   .. code-block:: none

      $ git add <some-path/file-name>

#. Commit changes locally with either of the following commands depending on the
   detail needed in the commit message. See :ref:`CommitMessages` for guidelines
   when writing longer commit messages.

   a. For short commit messages:

      .. code-block:: none

         $ git commit -m "Short Description of the Commit"

   --- OR ---

   b. For longer commit messages:

      .. code-block:: none

         $ git commit

      This opens the default system editor or the editor set in ``.gitconfig``,
      see :ref:`GitConfig` for more information.

#. Repeat the above modify, stage, and commit steps as needed to complete the
   work.

#. During the development cycle it is a good practice to periodically push local
   commits to the remote repo. To push your locally committed changes use:

   .. code-block:: none

      $ git push

   If changes exist on the remote copy of the branch that are not in your local
   version, Git will not allow you to push. You must first pull in the remote
   changes as described below and then push your local changes.

#. If you are collaborating with other developers on a feature you may need to
   update your local copy of a branch with changes from the remote version. To
   pull in remote changes use:

   .. code-block:: none

      $ git pull

   This will attempt to merge any remote commits into your local branch. If
   conflicts arise, Git will tell you. Follow the instructions in
   :ref:`ResolvingConflicts` to manually finish the merge. More experienced Git
   users may wish to use the ``--rebase`` option with ``git pull``. See
   :ref:`Rebasing` for more information.

   NOTE: Locally tracked changes must be committed (or stashed) prior to
   pulling remote changes.

#. As changes are merged into the PARENT branch, you may need to incorporate
   those changes into your working branch. To pull in changes from the remote
   PARENT branch use:

   .. code-block:: none

      $ git pull origin PARENT

   This will attempt to merge any changes on the remote PARENT branch into your
   local branch. If conflicts arise, Git will tell you. Follow the instructions
   in :ref:`ResolvingConflicts` to manually finish the merge. More experienced
   Git users may wish to use the ``--rebase`` option with ``git pull``. See
   :ref:`Rebasing` for more information.

   NOTE: Locally tracked changes must be committed (or stashed) prior to
   pulling remote changes.

#. Once the development work is complete make sure all changes have been pushed
   to the remote branch and open a pull request (PR) to start a code review for
   integrating the changes into the PARENT branch. See :ref:`PullRequests` for
   more details on the opening a PR and the review process.

.. _CommitMessages:

Commit Messages
^^^^^^^^^^^^^^^

The desired format for longer commit messages (more than a single line) is a
short descriptive title followed by a blank line, and then a detailed commit
message. For example, a commit making several changes to the ARKODE
initialization function might have the following message:

.. code-block:: none

   Retain history, remove LS/NLS init flag, skip init on reset

   * move initializations inside FIRST_INIT
   * retain error history on reset
   * clear error/step history in SetInitStep
   * skip stepper init checks on reset
   * remove updates to ark LS and NLS init functions
   * remove prototypes for old reinit utility function
   * update docs

.. _ResolvingConflicts:

Resolving Conflicts
^^^^^^^^^^^^^^^^^^^

When updating the local copy of a branch with changes from the remote version or
the PARENT branch merge conflicts may arise. In this case the conflicts can be
resolved using either of the following approaches.

Manually resolve the conflicts in an editor:

#. View the status of the branch for more info on the files with conflicts:

   .. code-block:: none

      $ git status

#. Use your editor of choice to resolve the conflicts

#. Add the resolved files:

   .. code-block:: none

      $ git add <path/file-name>

#. Repeat the above steps until all conflicts are resolved

#. If conducting a merge, commit the merge changes with

   .. code-block:: none

      $ git commit

   or if rebasing, continue the rebase with

   .. code-block:: none

      $ git rebase --continue

#. If desired, push the updated local branch to the remote repo with

   .. code-block:: none

      $ git push

Alternatively, if you have setup a merge tool in your ``.gitconfig`` file (see
:ref:`GitConfig` for more details), then use resolve the conflicts with the
merge tool:

#. Start the mergetool

   .. code-block:: none

      $ git mergetool

   Git will open each conflicted file in the merge tool one at a time. Resolve
   the conflicts, save the merged file, and close the merge tool. Git will then
   open the next conflicted file. As files are resolved and saved they will be
   automatically added to the merge/rebase commit.

#. When all conflicts are resolved, either commit the changes to the local
   branch if merging with

   .. code-block:: none

      $ git commit

   or if rebasing, continue the rebase with

   .. code-block:: none

      $ git rebase --continue

#. If desired, push the updated local branch to the remote repo with

   .. code-block:: none

      $ git push

.. _Rebasing:

Rebasing
^^^^^^^^

Rebasing is an alternative to merging that replays commits one at a time on top
of the HEAD of another branch (e.g., the remote copy of a local branch or the
branch's parent). This rewrites the commit history by reapplying old changes in
a new commit to keep the history linear and avoid a merge commit. Rewriting
history is potentially dangerous and should be done with caution.

As mentioned above when pulling changes the ``--rebase`` option can be used to
perform a rebase rather than a merge e.g.,

.. code-block:: none

   $ git pull --rebase

If conflicts arise during the rebase (Git will tell you), follow the steps in
:ref:`ResolvingConflicts` to resolve them. If you choose to continue working
locally and later run ``git pull --rebase`` again, you will likely have to
re-resolve conflicts that occurred in the previous rebase. In this case
``git rerere`` can be used to have Git remember how conflicts were resolved
previously.

Generally, rebasing changes on the PARENT branch e.g., when running
``git pull --rebase origin PARENT`` is not recommended as *all* of the commits
(local and remote) on the working branch will be replayed on the PARENT branch.
As such this will rewrite commits that have been pushed to the remote repo and
will therefore change the commit history that other developers have pulled and
are working on.
