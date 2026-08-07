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

.. _Contributing:

Contributing
============

There are two primary ways of contributing to SUNDIALS. The first way is by participating
in the development of SUNDIALS directly through contributions of code to the primary
`SUNDIALS repository <https://github.com/LLNL/sundials>`_. This is the best way to contribute
bug fixes and minor improvements.  At this time, the SUNDIALS team does not have the resources
to review and take in large additions to the code or significant new features.
Larger additions can be contributed as a SUNDIALS "addon" which is a component that may be
optionally downloaded by users and then compiled and installed with SUNDIALS.

.. note::

   New SUNDIALS developers should read the :ref:`Developer` guide for information on
   how to participate in SUNDIALS development.


Direct Contributions via Pull Requests
--------------------------------------

Direct contributions to SUNDIALS are made by opening a :ref:`Pull Request <OpenPR>`.

All new contributions to SUNDIALS must be made under the BSD 3-clause license.
The SUNDIALS team will not accept any file previously released under any other open
source license. By submitting code, the contributor gives irreversible consent
to the redistribution and/or modification of the contributed source code.

Please ensure that any pull request includes user guide additions or changes as
appropriate, has been tested, and includes a test for any added features.

Any added files submitted with a pull request must contain a header section at
the top including the originating author's name and file origin date as well as
a pointer to the SUNDIALS LICENSE and NOTICE files.

The act of submitting a pull request (with or without an explicit Signed-off-by
tag) will be understood as an affirmation of the following `Developer's
Certificate of Origin (DCO) <http://developercertificate.org/>`_.

.. code-block:: text

   Developer Certificate of Origin
   Version 1.1

   Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
   1 Letterman Drive
   Suite D4700
   San Francisco, CA, 94129

   Everyone is permitted to copy and distribute verbatim copies of this
   license document, but changing it is not allowed.


   Developer's Certificate of Origin 1.1

   By making a contribution to this project, I certify that:

   (a) The contribution was created in whole or in part by me and I
      have the right to submit it under the open source license
      indicated in the file; or

   (b) The contribution is based upon previous work that, to the best
      of my knowledge, is covered under an appropriate open source
      license and I have the right under that license to submit that
      work with modifications, whether created in whole or in part
      by me, under the same open source license (unless I am
      permitted to submit under a different license), as indicated
      in the file; or

   (c) The contribution was provided directly to me by some other
      person who certified (a), (b) or (c) and I have not modified
      it.

   (d) I understand and agree that this project and the contribution
      are public and that a record of the contribution (including all
      personal information I submit with it, including my sign-off) is
      maintained indefinitely and may be redistributed consistent with
      this project or the open source license(s) involved.


The DCO lets us know that you are entitled to contribute this code to SUNDIALS and that you are
willing to have it used in distributions and derivative works.

By including the DCO signature, you are stating that one or
more of the following is true of your contribution:

1.  You created this contribution/change and have the right to submit it
    to SUNDIALS; or
2.  You created this contribution/change based on a previous work with a
    compatible open source license; or
3.  This contribution/change has been provided to you by someone who did
    1 or 2 and you are submitting the contribution unchanged.
4.  You understand this contribution is public and may be redistributed as
    open source software under the BSD license.

All commits submitted to the SUNDIALS project need to have the following sign
off line in the commit message:

.. code-block:: text

   Signed-off-by: Jane Doe <jdoe@address.com>


Replacing Jane Doe’s details with your name and email address.

If you've set ``user.name`` and ``user.email`` in your Git configuration, you can
automatically add a sign off line at the end of the commit message by using the
``-s`` option (e.g., ``git commit -s``).


Contributions via SUNDIALS Addons
---------------------------------

SUNDIALS "addons" are community developed code additions for SUNDIALS that can be subsumed by the
SUNDIALS build system so that they have full access to all internal SUNDIALS symbols.
The intent is for SUNDIALS addons to function as if they are part of the SUNDIALS library,
while allowing them to potentially have different licenses
(although we encourage BSD-3-Clause still), code style
(although we encourage them to follow the SUNDIALS style outlined :ref:`here <SourceCode>`),
and they **are not maintained by the SUNDIALS team**.

Creating an addon
^^^^^^^^^^^^^^^^^

To create a SUNDIALS addon and use it there are a few things you need to do:

1. In your addon project, ensure that you have a ``CMakeLists.txt`` that uses the
   ``sundials_add_library`` CMake macro to create the library target. The best thing to do is simply
   copy from, or refer to, a ``CMakeLists.txt`` in the SUNDIALS ``src`` directory.
2. Follow the steps in the ``README.md`` file in the ``external/`` directory in the root of the SUNDIALS
   source code.

An example addon is available `here <https://github.com/sundials-codes/sundials-addon-example>`_.

Friends of SUNDIALS
^^^^^^^^^^^^^^^^^^^

The SUNDIALS team maintains a list of some SUNDIALS addons in our `Friends of SUNDIALS
<https://github.com/sundials-codes/friends-of-sundials>`_ repository. These addons are not
maintained by the SUNDIALS team, but have been developed in consultation with us.
**Currently we are only adding projects for existing collaborations**.
Please contact the development team if you are interested in collaborating on an addon.
