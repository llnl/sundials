..
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
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
   ----------------------------------------------------------------

.. _Publications:

#####################
SUNDIALS Publications
#####################

General
=======


High-level publications describing SUNDIALS in general include :cite:p:`RARBGW:26, GRWB:22, HBGLSSW:05`:

.. code-block:: latex

   @article{roberts2026sundials,
     title     = {New {Time} {Integrators} and {Capabilities} in {SUNDIALS} {Versions} 6.2.0-7.4.0},
     author    = {Roberts, Steven B and Aggul, Mustafa and Reynolds, Daniel R and Balos, Cody J and Gardner, David J and Woodward, Carol S},
     journal   = {ACM Transactions on Mathematical Software (TOMS)},
     publisher = {ACM},
     volume    = {52},
     number    = {1},
     pages     = {8:1--8:14},
     year      = {2026},
     doi       = {10.1145/3797888}
   }

.. code-block:: latex

   @article{gardner2022sundials,
     title     = {Enabling new flexibility in the {SUNDIALS} suite of nonlinear and differential/algebraic equation solvers},
     author    = {Gardner, David J and Reynolds, Daniel R and Woodward, Carol S and Balos, Cody J},
     journal   = {ACM Transactions on Mathematical Software (TOMS)},
     publisher = {ACM},
     volume    = {48},
     number    = {3},
     pages     = {1--24},
     year      = {2022},
     doi       = {10.1145/3539801}
   }

.. code-block:: latex

   @article{hindmarsh2005sundials,
     title     = {{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
     author    = {Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
     journal   = {ACM Transactions on Mathematical Software (TOMS)},
     publisher = {ACM},
     volume    = {31},
     number    = {3},
     pages     = {363--396},
     year      = {2005},
     doi       = {10.1145/1089014.1089020}
   }

GPU features in SUNDIALS are described in :cite:p:`balos2021enabling`:

.. code-block:: latex

   @article{balos2021enabling,
     title     = {{Enabling GPU accelerated computing in the SUNDIALS time integration library}},
     author    = {Balos, Cody J and Gardner, David J and Woodward, Carol S and Reynolds, Daniel R},
     journal   = {Parallel Computing},
     publisher = {Elsevier},
     volume    = {108},
     pages     = {102836},
     year      = {2021},
     doi       = {10.1016/j.parco.2021.102836}
   }


ARKODE
======

While ARKODE was included in :cite:p:`GRWB:22`, it was officially introduced in :cite:p:`RGWC:23`:

.. code-block:: latex

   @article{reynolds2023arkode,
     title     = {{ARKODE}: {A} flexible {IVP} solver infrastructure for one-step methods},
     author    = {Reynolds, Daniel R and Gardner, David J and Woodward, Carol S and Chinomona, Rujeko},
     journal   = {ACM Transactions on Mathematical Software (TOMS)},
     publisher = {ACM},
     volume    = {49},
     number    = {2},
     year      = {2023},
     pages     = {1--26},
     doi       = {10.1145/3594632}
   }

Time adaptivity for ARKODE's multirate solvers is described in :cite:p:`RAML:26`:

.. code-block:: latex

   @article{reynolds2026efficient,
     title     = {Efficient and {Flexible} {Multirate} {Temporal} {Adaptivity}},
     author    = {Reynolds, Daniel R. and Amihere, Sylvia and Mitchell, Dashon and Luan, Vu Thai},
     journal   = {Journal of Computational and Applied Mathematics},
     publisher = {Elsevier},
     year      = {2026},
     pages     = {117773},
     doi       = {10.1016/j.cam.2026.117773}
   }
