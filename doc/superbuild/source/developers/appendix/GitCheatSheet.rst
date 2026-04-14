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

.. _GitCheatSheet:

Git Cheat Sheet
===============

This section contains details helpful Git commands. For information on setting
up Git see the :ref:`GitSetup` section and for the SUNDIALS Git workflow see
the :ref:`Workflow` section.

**Command Help**

  Print a list of common Git commands used in various situations

  .. code-block:: none

     $ git --help

  Print full man page for a Git command

  .. code-block:: none

     $ git <git-command> --help

  Print list of available options for a Git command

  .. code-block:: none

     $ git <git-command> -h

**Graphical Tools**

  Availability depends the Git installation

  Visualize repo history

  .. code-block:: none

     $ gitk

  Graphical interface to Git

  .. code-block:: none

     $ git gui

**Branch**

  List all local branches (* denotes the current branch)

  .. code-block:: none

     $ git branch

  List all local branches with additional information

  .. code-block:: none

     $ git branch -vv

  List all remote branches

  .. code-block:: none

     $ git branch -r

  Create a new branch from an existing reference point (e.g., branch, commit
  hash, etc.)

  .. code-block:: none

     $ git branch <new-branch-name> <reference>

  Create a new branch from an existing reference point (e.g., branch, commit
  hash, etc.) and immediately checkout the new branch

  .. code-block:: none

     $ git checkout -b <new-branch-name> <reference>

**Add**

  Stage local changes to a file

  .. code-block:: none

     $ git add <file-name>

  Stage changes from a file in chunks

  .. code-block:: none

     $ git add -p <file-name>

  Unstage a committed file

  .. code-block:: none

     $ git reset HEAD -- <file-name>

  Discard local changes to a file and get the previously committed version

  .. code-block:: none

     $ git checkout -- <file-name>


  Discard local changes to a file in chunks

  .. code-block:: none

     $ git checkout -p <file-name>


**Commit**

  Commit staged files, this opens an editor for writing a commit message, an
  empty commit message will abort the commit

  .. code-block:: none

     $ git commit

  The desired format for a commit message is a short descriptive title
  followed by a blank line, and then a detailed commit message. For
  example, a commit making several changes to the ARKODE initialization
  function might have the following message:

  .. code-block:: none

     Retain history, remove LS/NLS init flag, skip init on reset

     * move initializations inside FIRST_INIT
     * retain error history on reset
     * clear error/step history in SetInitStep
     * skip stepper init checks on reset
     * remove updates to ark LS and NLS init functions
     * remove prototypes for old reinit utility function
     * update docs

  Commit staged files with short message

  .. code-block:: none

     $ git commit -m "short commit message"

  Amend the most recent commit (assuming it has not been pushed) to include the
  changes in this commit

  .. code-block:: none

     $ git commit --amend

  This is useful for adding forgotten changes to the last commit or for editing
  the commit message in the last commit (if no files are staged).

**Push and Pull**

  Push a new branch and its commits to the remote repository and set the
  upstream tracking branch

  .. code-block:: none

     $ git push -u origin <branch-name>

  Push committed changes on the current branch to the remote repository

  .. code-block:: none

     $ git push

  Note: This is the same as ``git push origin <current-branch-name>``.

  Fetch and merge remote changes for the current branch

  .. code-block:: none

     $ git pull

  Note: This is the same as ``git pull origin <current-branch-name>``.

  Fetch and merge changes from a specific remote branch into the current branch

  .. code-block:: none

     $ git pull origin <branch-name>

  Fetch and rebase remote changes for the current branch

  .. code-block:: none

     $ git pull --rebase

  Note: Use caution when rebasing, see :ref:`Rebasing` for more details.

  Fetch and rebase changes from a specific remote branch into the current branch

  .. code-block:: none

     $ git pull --rebase origin <branch-name>

  Note: Use caution when rebasing, see :ref:`Rebasing` for more details.

**Merge**

  Merge a different local branch to your current branch

  .. code-block:: none

     $ git merge <branch-name>

  Resolve merge conflicts with a visual diff/merge tool

  .. code-block:: none

     $ git mergetool

  Find the newest common ancestor (fork point) of two reference points
  (branches, commits, etc.)

  .. code-block:: none

     $ git merge-base <reference1> <reference2>

**Rebase**

  Interactively rebase the committed but not pushed changes on the current
  branch

  .. code-block:: none

     $ git rebase -i

  Note: This command is useful for cleaning up the local commits history (e.g.,
  reordering, squashing, updating commit messages, etc.) before pushing to the
  remote repository.

  Rebase the committed but not pushed changes on the current branch and execute
  a given command (``cmd``) between each step in the rebase

  .. code-block:: none

     $ git rebase -x "cmd"

  Note: This command can be useful for debugging commits e.g., ``cmd`` can be
  used to run the test suite after each commit to find which commit causes the
  test suite to fail.

**Cherry Pick**

  Apply a specific commit to the current branch

  .. code-block:: none

     $ git cherry-pick <commit-hash>

**Status and Differences**

  Print information on current local repo status including unstaged (changed and
  not added) files, staged (changed and added) files, and untracked files.

  .. code-block:: none

     $ git status

  Show *all* differences between unstaged changes and the current HEAD

  .. code-block:: none

     $ git diff

  Show the differences between unstaged changes and the current HEAD for a
  specific file

  .. code-block:: none

     $ git diff <file-name>

  Show differences between *all* staged files and current HEAD

  .. code-block:: none

     $ git diff --staged

  List *all* files changed in the current branch compared to different reference
  point (branch, commit hash, etc.)

  .. code-block:: none

     $ git diff --name-only <reference>

  Compare files between two branches

  .. code-block:: none

     $ git diff <branch1>..<branch2> -- <file-name>

  To view the differences going from the remote file to the local file

  .. code-block:: none

     $ git diff remotename/branchname:remote/path/file1 local/path/file1

  To view the differences going from the local file to the remote file

  .. code-block:: none

     $ git diff HEAD:local/path/file1 remotename/branchname:remote/path/file1

  To view the differences between files at any two reference points (e.g.
  branches, commit hashes, etc.)

  .. code-block:: none

     $ git diff ref1:path/to/file1 ref2:path/to/file2

  Note: In the above commands ``diff`` can be replaced with ``difftool`` if a
  visual diff tool has been setup with Git.

**Log**

  Show commit log

  .. code-block:: none

     $ git log

  Show commit log for the last n commits

  .. code-block:: none

     $ git log -n

  Show commit log with more change information

  .. code-block:: none

     $ git log --stat

  Show all commits impacting a specific file

  .. code-block:: none

     $ git log <file-name>

**Stash**

  Save uncommitted changes in the stash

  .. code-block:: none

     $ git stash

  Save uncommitted changes in the stash with a stash message

  .. code-block:: none

     $ git stash save "stash message"

  View saved changes in the stash

  .. code-block:: none

     $ git stash list

  Apply to first set of stashed changes and remove the changes from the stash

  .. code-block:: none

     $ git stash pop

  Apply changes from the stash (does not remove changes from the stash)

  .. code-block:: none

     $ git stash apply <stash-name>

  Remove changes from the stash

  .. code-block:: none

     $ git stash drop <stash-name>

  Show difference between current HEAD and stashed work

  .. code-block:: none

     $ git stash show -p <stash-name>
