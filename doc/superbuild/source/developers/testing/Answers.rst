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

.. _Answers:

Answer Files
============

When solving ODEs/DAEs and nonlinear algebraic equations numerically, there is not always a known
analytical solution for the problem to compare with. Additionally, there may be many paths to a
correct solution and some may be less efficient than others. This presents an issue for testing
SUNDIALS. How do we ensure SUNDIALS is working correctly if we cannot write down what "correct" is?
To solve this problem SUNDIALS has tests (which often also serve as example codes) that have been
verified by the author to be functioning correctly. These tests produce output which is stored in a
file known as the "answer" file (aka "out" or "output" file). The output typically consist of the
solution and/or some quantities derived from the solution as well as statistics about how the
integration and/or solve (e.g., number of time steps). When changes are made to SUNDIALS, we ensure
that these output files do not change unless it is expected/desired. Practically, this is ensured
by our CI.

Due to difference in microarichtectures and the nature of floating point arithmetic, it is quite
possible that the output generated on one machine may differ from the output generated on another.
As such, we specify that the the answer files that are embedded in ``examples/`` (the ``.out``
files) should match what is produced on the Jenkins CI machine.  We also have a `separate git
repository <https://github.com/sundials-codes/answers>`_ which holds answers for other  machines,
such as the GitHub Actions virtual machines.


Updating and Adding Answer Files for Existing Machines
------------------------------------------------------

Sometimes it is necessary to update the answer files, or add a new file, so that the tests pass, if a change in output
is expected/desired. Changing output files requires careful verification that the output is still
"correct". The procedure for updating answer files is as follows:

- For the ``.out`` files embedded in examples, simply overwrite the file with the new output. If the
  tests still fail, then you will need to download the output file generated on the Jenkins CI
  machine and commit it. Extracting this file from the Jenkins machine may require help from David
  or Cody currently.

- For the ``.out`` files in `sundials-codes/answers <https://github.com/sundials-codes/answers>`_:

  #. Navigate to `test/answers` in the sundials repository. If there is nothing in this directory,
     you need to run ``git submodule init && git submodule update``. Alternatively clone the answers
     repository outside of the SUNDIALS source tree.

  #. Checkout the main branch: ``git checkout main``.

  #. Create a new branch off of main. It is ideal if you use a branch name that matches the name of
     the branch that you are developing in sundials on: ``git checkout -b <branch name>``.

  #. Update the relevant ``.out`` files in
     ``test/answers/linux-ubuntu20.04-x86_64/gcc-9.4.0/<double|single|extended>`` (this is the files
     used for the GitHub Actions CI). Like with the embedded ``.out`` files, you can try and use the
     output generated on your machine, but it may be different enough from what is produced on the
     GitHub actions to trigger a failure. You can download the output files generated on the GitHub
     machines by going to `<https://github.com/LLNL/sundials/actions>`_ finding your failing
     test, clicking it, then at the bottom downloading the "output_files" artifact (you can also find your
     failing test at the bottom of a PR). The downloaded zip file will be the SUNDIALS build
     directory. To update just the failed files you can use ``scripts/updateOutFiles.py`` e.g.,

     .. code:: console

        cd scripts
        ./updateOutFiles.py <source path> <destination path>

     When updating files in the answers repo the source path is the path to the build directory
     downloaded from GitHub and the destination path is ``<path to answers
     repo>/linux-ubuntu20.04-x86_64/gcc-9.4.0/<double|single|extended>``. When updating with output
     from Jenkins the source path is the build directory for the Jenkins run and the destination
     path is the top level of the sundials repo (i.e., ../. if running from scripts). The script
     will automatically traverse the examples and test/unit_test directories to update output files
     for failed tests.

  #. If you need to use ssh authentication for pushing to GitHub, then you may need to update the
     remote for the submodule:

   .. code:: console

      git remote set-url origin git@github.com:sundials-codes/answers.git

  #. Create a pull-request in `sundials-codes/answers <https://github.com/sundials-codes/answers>`_
     with your updates.


Adding Answer Files for a New Machine
-------------------------------------

- For the ``.out`` files embedded in examples, simply overwrite the file with the new output. If the
  tests running on Jenkins still fail, then you will need to download the output file generated on
  the Jenkins CI machine and commit it. Extracting this file from the Jenkins machine may require
  help from David or Cody currently.

- For the ``.out`` files in `sundials-codes/answers <https://github.com/sundials-codes/answers>`_:

  #. To add new answers, start by copying the relevant answers in ``linux-ubuntu20.04-x86_64`` to
     your new directory (following the naming scheme outlined above).

  #. Next Commit these files and open a PR to the "staging" branch.

  #. Once the PR is merged, you can generate your new answer files and overwrite the ones you copied
     from ``linux-ubuntu20.04-x86_64``.

  #. Compare the diff and make sure it is reasonable, then commit.

  #. Finally, open a new PR targeting the "staging" branch.

  #. Eventually, "staging" will be merged into main.
