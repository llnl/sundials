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

.. _GitSetup:

Git Setup
=========

.. _SSHKeys:

SSH Keys
--------

One-time setup of ssh keys to avoid needing a password with each repository
transaction.

#. Create your ssh key using the ssh-keygen utility:

   a. Skip this first command if you already have a public/private ssh key
      (i.e., you already have an ``id_rsa.pub`` file in your ``~/.ssh``
      directory):

      .. code-block:: none

        $ ssh-keygen -t rsa

   b. Make sure that your ssh directory (``~/.ssh``) is readable only by you.

      .. code-block:: none

        $ chmod -R go-rwx ~/.ssh

#. Install the key in GitHub following the instructions `here <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_

.. _GitConfig:

Git Configuration
-----------------

Global Git configuration information is stored in ``~/.gitconfig`` and can be
set with the ``git config --global`` command or by editing the ``.gitconfig``
file manually. More information is available in the `git config documentation
<https://git-scm.com/docs/git-config>`_. Note in the comments in the
``.gitconfig`` file begin with ``#``.

#. Tell Git who you are:

   .. code-block:: none

     $ git config --global user.name "Your Name"
     $ git config --global user.email "yourname@example.com"

   --- OR ---

   .. code-block:: none

     [user]
         name = Your Name             # name recorded in commits
         email = yourname@example.com # email recorded in commits

#. Tell Git to only push changes on your current branch and only if the upstream
   branch has the same name i.e., don't default to push changes to all branches
   or branches with a different name. Note the behavior is the default beginning
   with Git v2.0.

   .. code-block:: none

      $ git config --global push.default simple

   --- OR ---

   .. code-block:: none

      [push]
          default = "simple" # push current branch to upstream branch with same name

#. Tell Git about the .git-blame-ignore-revs file.

   git config blame.ignoreRevsFile .git-blame-ignore-revs

#. Tell Git which editor you want to use by default (e.g., Emacs):

   .. code-block:: none

      $ git config --global core.editor "emacs"

   --- OR ---

   .. code-block:: none

      [core]
           editor = emacs # default editor

#. Enable color output in Git

   .. code-block:: none

      $ git config --global color.ui "true"

   --- OR ---

   .. code-block:: none

      [color]
           ui = true # enable color output

The following settings enable using a graphical diff tool to resolve conflicts
during a merge or when viewing diffs between files. These settings are optional
but may be useful. The settings below are given for the meld diff tool. Similar
settings can be used with emerge, gvimdiff, kdiff3, vimdiff, and tortoisemerge.

#. To add a merge tool (invoked by ``git mergetool``), add the following
   to the ``~/.gitconfig`` file:

   .. code-block:: none

      [mergetool "meld"]
          # command to invoke the merge tool for newer versions of
          # meld which use the '--output' option
          cmd = meld "$LOCAL" "$MERGED" "$REMOTE" --output "$MERGED"

#. To add a diff tool (invoked by ``git difftool``), add the following to the
   ``~/.gitconfig`` file:

   .. code-block:: none

      [diff]
          # which diff tool Git should use
          tool = meld
      [difftool]
          # do not prompt before each invocation of the diff tool
          prompt = false
      [difftool "meld"]
          # command to invoke the diff tool
          cmd = meld "$LOCAL" "$REMOTE"

Additionally, Git provides helpful scripts to enable auto-completion of Git
commands and to display the current status in the command line prompt for
various shells. The scripts ``git-completion.*`` and ``git-prompt.sh`` can be
obtained from the `contrib/completion directory
<https://github.com/git/git/tree/master/contrib/completion>`_ in the Git source
repository on GitHub.

For example with Bash, auto-completion can be enabled by adding

.. code-block:: bash

   source <some-path>/git-completion.bash

to your ``.bashrc`` file. Similarly displaying the Git command line prompt
information can be enabled by adding

.. code-block:: bash

   export GIT_PS1_SHOWDIRTYSTATE="true"       # unstaged *, staged +
   export GIT_PS1_SHOWSTASHSTATE="true"       # stashed $
   export GIT_PS1_SHOWUNTRACKEDFILES="true"   # untracked %
   export GIT_PS1_SHOWUPSTREAM="auto verbose" # ahead +, behind -, diverged +-, same =
   export GIT_PS1_SHOWCOLORHINTS="true"
   source <some-path>/git-prompt.sh

   export PROMPT_COMMAND='__git_ps1 "[$(date +%k:%M:%S)] \u@\h \w" "\n$"'

to your ``.bashrc`` file.

.. _CloneRepo:

Cloning the Repository
----------------------

To clone a copy of the SUNDIALS repository use one of the following commands:

#. Clone the repository with SSH Keys:

   .. code-block:: none

      $ git clone --recurse-submodules git@github.com:LLNL/sundials.git

--- OR ---

#. Clone the repository with https (requires authenticating with your username and password or access token to push)

   .. code-block:: none

      $ git clone --recurse-submodules https://github.com/LLNL/sundials.git

After cloning the repository you will be on the ``main`` branch by default
however, the ``develop`` and ``main`` branches are protected branches and can
not be updated directly. In order to make changes to either of these branch you
must create a new branch from ``develop``, make the desired modifications, and
issue a pull request to have the changes merged into the parent branch. See the
:ref:`Workflow` section for more information.
