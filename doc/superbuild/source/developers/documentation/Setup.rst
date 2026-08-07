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

.. _Documentation.Setup:

Getting Started
===============

To build the documentation with Sphinx you will need `Python
<https://www.python.org/>`__ version 3.9 or newer. Generating the PDF versions
of the documentation requires a `LaTeX <https://www.latex-project.org>`__
installation and building the developer documentation requires `Graphviz
<https://graphviz.org/>`__ for generating flowcharts. All of these are generally
available from Linux package managers or `Homebrew <https://brew.sh/>`__ on Mac.

Setting Up Sphinx
-----------------

Once you have the dependencies installed, Sphinx and the necessary extensions
can be installed using the requirements file, ``doc/requirements.txt``, and
running the pip command below from the ``doc`` directory.

.. code-block:: shell

   pip install -r requirements.txt

A Python virtual environment can also be used by activating the environment
before running the above commands.

Building the Documentation
--------------------------

The HTML or PDF versions the package user guides can be built individually or an
HTML "super build" of the documentation can be generated that includes all of
the package guides along with the developer guide.

To generate the HTML files for the super build, run the following:

.. code-block:: shell

   cd doc/superbuild
   make html

To view the documentation, open the ``build/html/index.html`` file in a web
browser e.g., using ``start`` (Windows), ``open`` on Max, or ``xdg-open`` on
Linux.

If you want to build the HTML or PDF user guide for an individual package, run
the following:

.. code-block:: shell

   cd doc/<package>/guide
   make <html|latexpdf>

The example guide for ARKODE can be built with

.. code-block:: shell

   cd doc/arkode/examples
   make <html|latexpdf>

and the other package example guides can be built with

.. code-block:: shell

   cd doc/<package>
   make ex

Organization
------------

As noted above, each package in SUNDIALS has a directory within the ``doc``
folder and the user and example guides can be compiled separately for each
package. Documentation shared across packages (e.g., vectors, matrices, solvers,
etc.) is located in the ``doc/shared`` directory. These files are included in
the package documentation source using the Sphinx ``.. include::``
directive. Generally speaking, symlinks to files are not compatible with Sphinx
and is why we utilize ``_link.rst`` files with the ``.. include::`` statements.
Common figures and other assets go under ``doc/shared/figs`` while package
specific figures go in the ``doc/shared/figs/<package>`` directory. All
citations go into ``doc/shared/sundials.bib``.
