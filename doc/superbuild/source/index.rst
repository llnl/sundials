..
   Programmer(s): Cody J. Balos @ LLNL
   ----------------------------------------------------------------
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
   ----------------------------------------------------------------

########
SUNDIALS
########

The `SUNDIALS <https://computing.llnl.gov/projects/sundials>`_ library of time
integrators and nonlinear solvers provides robust and efficient numerical
methods for ordinary differential equations (ODEs), differential-algebraic
equations (DAEs), and nonlinear algebraic systems. SUNDIALS is freely available
and developed on `GitHub <https://github.com/LLNL/sundials>`_.

SUNDIALS is comprised of following packages:

* :ref:`ARKODE <ARKODE>`, a solver with one-step methods for stiff, nonstiff,
  mixed stiff-nonstiff, and multirate ODE systems.

* :ref:`CVODE <CVODE>`, a solver with Adams and BDF methods for stiff and
  nonstiff ODE systems.

* :ref:`CVODES <CVODES>`, an extension of CVODE with forward and adjoint
  sensitivity analysis capabilities for stiff and nonstiff ODE systems.

* :ref:`IDA <IDA>`, a solver with BDF methods for DAE systems.

* :ref:`IDAS <IDAS>`, an extension of IDA with forward and adjoint sensitivity
  analysis capabilities for DAE systems.

* :ref:`KINSOL <KINSOL>`, a solver for nonlinear algebraic systems.

The SUNDIALS packages share many components and are organized as a family built
on a common infrastructure including abstract interfaces for vectors, matrices,
and algebraic solvers. Several implementations of these interfaces are provided
with SUNDIALS supporting a range of parallel computing paradigms including
shared-memory, distributed memory, and GPU computing.

Citing
======

.. include:: ../../shared/cite_sundials.rst

When using the ARKODE package from SUNDIALS, please also cite:

.. code-block:: latex

   @article{reynolds2023arkode,
     title   = {{ARKODE: A flexible IVP solver infrastructure for one-step methods}},
     author  = {Reynolds, Daniel R and Gardner, David J and Woodward, Carol S and Chinomona, Rujeko},
     journal = {ACM Transactions on Mathematical Software},
     volume  = {49},
     number  = {2},
     pages   = {1--26},
     year    = {2023},
     doi     = {10.1145/3594632}
   }

We also ask that users cite the documentation for the package and version that
they are using rather than the combined SUNDIALS online guide:

.. parsed-literal::

   @Misc{arkodeDocumentation,
      author       = {Daniel R. Reynolds and David J. Gardner and Carol S. Woodward, Rujeko Chinomona and Cody J. Balos},
      title        = {User Documentation for ARKODE},
      year         = {|YEAR|},
      note         = {|ARKODE_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/arkode},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/arkode}}
   }

.. parsed-literal::

   @Misc{cvodeDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for CVODE},
      year         = {|YEAR|},
      note         = {|CVODE_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/cvode},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/cvode}}
   }

.. parsed-literal::

   @Misc{cvodesDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for CVODES},
      year         = {|YEAR|},
      note         = {|CVODES_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/cvodes},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/cvodes}}
   }

.. parsed-literal::

   @Misc{idaDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for IDA},
      year         = {|YEAR|},
      note         = {|IDA_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/ida},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/ida}}
   }

.. parsed-literal::

   @Misc{idasDocumentation,
      author       = {Radu Serban and Cosmin Petra and Alan C. Hindmarsh and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for IDAS},
      year         = {|YEAR|},
      note         = {|IDAS_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/idas},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/idas}}
   }

.. parsed-literal::

   @Misc{kinsolDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for KINSOL},
      year         = {|YEAR|},
      note         = {|KINSOL_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/kinsol},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/kinsol}}
   }

Contributors
============

The SUNDIALS library has been developed over many years by a number of
contributors. The current SUNDIALS team consists of Cody J.  Balos,
David J. Gardner, Alan C. Hindmarsh, Daniel R. Reynolds, and Carol S.
Woodward. We thank Radu Serban for significant and critical past contributions.

Other contributors to SUNDIALS include: Mustafa Aggul, James Almgren-Bell, Lawrence E. Banks,
Peter N. Brown, George Byrne, Rujeko Chinomona, Scott D. Cohen, Aaron Collier,
Keith E. Grant, Steven L. Lee, Shelby L. Lockhart, John Loffeld, Daniel McGreer,
Yu Pan, Slaven Peles, Cosmin Petra, Steven B. Roberts, H. Hunter Schwartz,
Jean M. Sexton, Dan Shumaker, Steve G. Smith, Shahbaj Sohal, Allan G. Taylor,
Hilari C. Tiedeman, Chris White, Ting Yan, and Ulrike M. Yang.

Acknowledgments
===============

This material is based on work supported by the U.S. Department of Energy,
Office of Science, Office of Advanced Scientific Computing Research, Scientific
Discovery through Advanced Computing (SciDAC) program via the Frameworks,
Algorithms, and Scalable Technologies for Mathematics (FASTMath) Institute under
DOE awards DE-AC52-07NA27344 and DE-SC-0021354.

This material is also based on work supported by the U.S. Department of Energy,
Office of Science, Office of Advanced Scientific Computing Research,
Next-Generation Scientific Software Technologies program under contract
DE-AC52-07NA27344.  Additional support is also provided by SciDAC
partnerships with the U.S. Department of Energy’s FES, NP, BES, OE, and BER
offices as well as the LLNL Institutional Scientific Capability Portfolio.

SUNDIALS License and Notices
============================

.. include:: ../../shared/LicenseReleaseNumbers.rst


.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :numbered:
   :hidden:

   sundials/index.rst
   arkode/index.rst
   cvode/index.rst
   cvodes/index.rst
   ida/index.rst
   idas/index.rst
   kinsol/index.rst
   nvectors/index.rst
   sunmatrix/index.rst
   sunlinsol/index.rst
   sunnonlinsol/index.rst
   sundomeigest/index.rst
   sunadaptcontroller/index.rst
   sunstepper/index.rst
   sunadjoint/index.rst
   sunmemory/index.rst
   History_link.rst
   Changelog_link.rst
   FAQ_link.rst

.. toctree::
   :caption: INTERFACES
   :maxdepth: 1
   :numbered:
   :hidden:

   Fortran/index.rst

.. toctree::
   :caption: EXAMPLES
   :maxdepth: 1
   :hidden:

   examples/arkode/index-superbuild.rst
   examples/cvode/index-superbuild.rst
   examples/cvodes/index-superbuild.rst
   examples/ida/index-superbuild.rst
   examples/idas/index-superbuild.rst
   examples/kinsol/index-superbuild.rst

.. toctree::
   :caption: DEVELOPMENT
   :maxdepth: 1
   :hidden:

   contributing/index.rst
   developers/index.rst

.. toctree::
   :caption: REFERENCES
   :maxdepth: 1
   :hidden:

   References
