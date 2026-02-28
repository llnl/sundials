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

.. _KINSOL.Examples.CXX:

C++ Example Problems
====================

Parallel Matrix-Free Example: kin_heat2D_nonlin_p
-------------------------------------------------

As an illustration using KINSOL for the solution of nonlinear systems in
parallel, we give a sample program called ``kin_heat2D_nonlin_p.cpp``. It uses
the KINSOL fixed-point (``KIN_FP``) iteration with Anderson Acceleration and the
:ref:`MPI parallel vector <NVectors.NVParallel>` for the solution of a
steady-state 2D heat equation with an additional nonlinear term defined by
:math:`c(u)`:

.. math::

   b = \nabla \cdot (D \nabla u) + c(u) \quad \text{in} \quad \mathcal{D} = [0,1] \times [0,1]

where :math:`D` is a diagonal matrix with entries :math:`k_x` and :math:`k_y`
for the diffusivity in the :math:`x` and :math:`y` directions, respectively. The
boundary conditions are

.. math::

   u(0,y) = u(1,y) = u(x,0) = u(x,1) = 0.

We chose the analytical solution to be

.. math::

   u_{\text{exact}} = u(x,y) = \sin^2(\pi x) \sin^2(\pi y)

Hence, we define the static term :math:`b` as follows

.. math::

   b = k_x 2 \pi^2 (\cos^2(\pi x) - \sin^2(\pi x)) \sin^2(\pi y)
       + k_y 2 \pi^2 (\cos^2(\pi y) - \sin^2(\pi y)) \sin^2(\pi x)
       + c(u_{\text{exact}})

The spatial derivatives are computed using second-order centered differences,
with the data distributed over :math:`n_x \times n_y` points on a uniform
spatial grid. The problem is set up to use spatial grid parameters
:math:`n_x=64` and :math:`n_y=64`, with heat conductivity parameters
:math:`k_x=1.0` and :math:`k_y=1.0`.

This problem is solved via a fixed point iteration with Anderson acceleration,
where the fixed point is function formed by adding :math:`u` to both sides of
the equation, i.e.,

.. math::

   b + u = \nabla \cdot (D \nabla u) + c(u) + u,

so that the fixed point function is

.. math::

   G(u) = \nabla \cdot (D \nabla u) + c(u) + u - b.

The problem is run using a tolerance of :math:`10^{-8}`, and a starting vector
containing all ones. This example highlights the use of the various
orthogonalization routine options within Anderson Acceleration, passed to the
example problem via the ``--orthaa`` flag. Available options include 0
(``KIN_ORTH_MGS``), 1 (``KIN_ORTH_ICWY``), 2 (``KIN_ORTH_CGS2``), and 3
(``KIN_ORTH_DCGS2``).

The following tables contain all available input parameters when running the
example problem.

Optional Input Parameter Flags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| Flag                          | Description                                                   |
+===============================+===============================================================+
| ``--mesh <nx> <ny>``          | mesh points in the x and y directions                         |
+-------------------------------+---------------------------------------------------------------+
| ``--np <npx> <npy>``          | number of MPI processes in the x and y directions             |
+-------------------------------+---------------------------------------------------------------+
| ``--domain <xu> <yu>``        | domain upper bound in the x and y direction                   |
+-------------------------------+---------------------------------------------------------------+
| ``--k <kx> <ky>``             | diffusion coefficients                                        |
+-------------------------------+---------------------------------------------------------------+
| ``--rtol <rtol>``             | relative tolerance                                            |
+-------------------------------+---------------------------------------------------------------+
| ``--maa <maa>``               | number of previous residuals for Anderson Acceleration        |
+-------------------------------+---------------------------------------------------------------+
| ``--damping <damping>``       | damping parameter for Anderson Acceleration                   |
+-------------------------------+---------------------------------------------------------------+
| ``--orthaa <orthaa>``         | orthogonalization routine used in Anderson Acceleration       |
+-------------------------------+---------------------------------------------------------------+
| ``--maxits <maxits>``         | max number of iterations                                      |
+-------------------------------+---------------------------------------------------------------+
| ``--c <cu>``                  | nonlinear function choice (integer between 1 - 17)            |
+-------------------------------+---------------------------------------------------------------+
| ``--timing``                  | print timing data                                             |
+-------------------------------+---------------------------------------------------------------+
| ``--help``                    | print available input parameters and exit                     |
+-------------------------------+---------------------------------------------------------------+

Input Parameter Flags for Nonlinear Function :math:`c(u)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+------------+--------------------------------------------------------------------------------------+
| Flag       | Function                                                                             |
+============+======================================================================================+
| ``--c 1``  | :math:`c(u) = u`                                                                     |
+------------+--------------------------------------------------------------------------------------+
| ``--c 2``  | :math:`c(u) = u^3 - u`                                                               |
+------------+--------------------------------------------------------------------------------------+
| ``--c 3``  | :math:`c(u) = u - u^2`                                                               |
+------------+--------------------------------------------------------------------------------------+
| ``--c 4``  | :math:`c(u) = e^u`                                                                   |
+------------+--------------------------------------------------------------------------------------+
| ``--c 5``  | :math:`c(u) = u^4`                                                                   |
+------------+--------------------------------------------------------------------------------------+
| ``--c 6``  | :math:`c(u) = \cos^2(u) - \sin^2(u)`                                                 |
+------------+--------------------------------------------------------------------------------------+
| ``--c 7``  | :math:`c(u) = \cos^2(u) - \sin^2(u) - e^u`                                           |
+------------+--------------------------------------------------------------------------------------+
| ``--c 8``  | :math:`c(u) = e^u u^4 - u e^{\cos(u)}`                                               |
+------------+--------------------------------------------------------------------------------------+
| ``--c 9``  | :math:`c(u) = e^{(\cos^2(u))}`                                                       |
+------------+--------------------------------------------------------------------------------------+
| ``--c 10`` | :math:`c(u) = 10(u - u^2)`                                                           |
+------------+--------------------------------------------------------------------------------------+
| ``--c 11`` | :math:`c(u) = -13 + u + ((5-u)u - 2)u`                                               |
+------------+--------------------------------------------------------------------------------------+
| ``--c 12`` | :math:`c(u) = \sqrt{5}(u - u^2)`                                                     |
+------------+--------------------------------------------------------------------------------------+
| ``--c 13`` | :math:`c(u) = (u - e^u)^2 + (u + u \sin(u) - \cos(u))^2`                             |
+------------+--------------------------------------------------------------------------------------+
| ``--c 14`` | :math:`c(u) = u + u e^u + u e^{-u}`                                                  |
+------------+--------------------------------------------------------------------------------------+
| ``--c 15`` | :math:`c(u) = u + u e^u + u e^{-u} + (u - e^u)^2`                                    |
+------------+--------------------------------------------------------------------------------------+
| ``--c 16`` | :math:`c(u) = u + u e^u + u e^{-u} + (u - e^u)^2 + (u + u\sin(u) - \cos(u))^2`       |
+------------+--------------------------------------------------------------------------------------+
| ``--c 17`` | :math:`c(u) = u + u e^{-u} + e^u (u + \sin(u) - \cos(u))^3`                          |
+------------+--------------------------------------------------------------------------------------+

Parallel Example Using hypre: kin_heat2D_nonlin_hypre_pfmg
----------------------------------------------------------

As an illustration of the use of the KINSOL package for the solution of
nonlinear systems in parallel and using *hypre* linear solvers, we give a sample
program called ``kin_heat2D_nonlin_hypre_pfmg.cpp``. It uses the KINSOL
fixed-point (``KIN_FP``) iteration with Anderson Acceleration and the :ref:`MPI
parallel vector <NVectors.NVParallel>` for the solution of a steady-state 2D
heat equation with an additional nonlinear term defined by :math:`c(u)`:

.. math::

   b = \nabla \cdot (D \nabla u) + c(u) \quad \text{in} \quad \mathcal{D} = [0,1] \times [0,1]

where :math:`D` is a diagonal matrix with entries :math:`k_x` and :math:`k_y`
for the diffusivity in the :math:`x` and :math:`y` directions, respectively. The
boundary conditions are

.. math::

   u(0,y) = u(1,y) = u(x,0) = u(x,1) = 0.

We chose the analytical solution to be

.. math::

   u_{\text{exact}} = u(x,y) = \sin^2(\pi x) \sin^2(\pi y)

Hence, we define the static term :math:`b` as follows

.. math::

   b = k_x 2 \pi^2 (\cos^2(\pi x) - \sin^2(\pi x)) \sin^2(\pi y)
       + k_y 2 \pi^2 (\cos^2(\pi y) - \sin^2(\pi y)) \sin^2(\pi x)
       + c(u_{\text{exact}})

The spatial derivatives are computed using second-order centered differences,
with the data distributed over :math:`n_x \times n_y` points on a uniform
spatial grid. The problem is set up to use spatial grid parameters
:math:`n_x=64` and :math:`n_y=64`, with heat conductivity parameters
:math:`k_x=1.0` and :math:`k_y=1.0`.

This problem is solved via a fixed point iteration with Anderson acceleration,
where the fixed point function is formed by implementing the Laplacian as a
matrix-vector product,

.. math::

   b = A u + c(u)

and solving for :math:`u` results in the fixed point function

.. math::

   G(u) = A^{-1} (b - c(u)).

The problem is run using a tolerance of :math:`10^{-8}`, and a starting vector
containing all ones. The linear system solve is executed using the :ref:`PCG
linear solver <SUNLinSol.PCG>` with the *hypre* PFMG preconditioner. The setup
of the linear solver can be found in the ``Setup_LS`` function, and setup of the
*hypre* preconditioner can be found in the ``Setup_Hypre`` function. This
example highlights the use of the various orthogonalization routine options
within Anderson Acceleration, passed to the example problem via the ``--orthaa``
flag. Available options include 0 (``KIN_ORTH_MGS``), 1 (``KIN_ORTH_ICWY``), 2
(``KIN_ORTH_CGS2``), and 3 (``KIN_ORTH_DCGS2``).

All input parameter flags available for the previous example are also available
for this problem. In addition, all runtime flags controlling the linear solver
and *hypre* related parameters are set using the flags in the following table.

Optional Input Parameter Flags for hypre
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+---------------------------------+---------------------------------------------------------------+
| Flag                            | Description                                                   |
+=================================+===============================================================+
| ``--lsinfo``                    | output residual history for PCG                               |
+---------------------------------+---------------------------------------------------------------+
| ``--liniters <liniters>``       | max number of iterations for PCG                              |
+---------------------------------+---------------------------------------------------------------+
| ``--epslin <epslin>``           | linear tolerance for PCG                                      |
+---------------------------------+---------------------------------------------------------------+
| ``--pfmg_relax <pfmg_relax>``   | relaxation type in PFMG                                       |
+---------------------------------+---------------------------------------------------------------+
| ``--pfmg_nrelax <pfmg_nrelax>`` | pre/post relaxation sweeps in PFMG                            |
+---------------------------------+---------------------------------------------------------------+
