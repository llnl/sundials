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

.. _Benchmarks.AdvectionReaction:


Advection-Reaction Benchmark
----------------------------

This benchmark problem implements a 3D advection-reaction equation using the
RAJA performance portability layer with serial, CUDA, or HIP backends.


Problem description
^^^^^^^^^^^^^^^^^^^

This code simulates the advection and reaction of three chemical species where
the reaction mechanism is a variation of the Brusselator problem from chemical
kinetics. The PDE system is given by

.. math::

   u_t &= -c \nabla u + A - (w+1) u + v u^2 \\
   v_t &= -c \nabla v + w u - v u^2 \\
   w_t &= -c \nabla w + (B - w) / \epsilon - w u

where :math:`u`, :math:`v`, and :math:`w` are chemical concentrations, :math:`c`
is the advection speed, :math:`A` and :math:`B` are the concentrations of
chemical species that remain constant over space and time, and :math:`\epsilon`
is a parameter that varies the stiffness of the system. The problem is solved on
the domain :math:`(x,y,z) \equiv \mathbf{x} \in [\mathbf{0}, \mathbf{x_{\text{max}}}]^3`,
for times :math:`t \in [0,t_f]`. The initial condition is

.. math::

   u(0,\mathbf{x}) &= A + p(\mathbf{x}) \\
   v(0,\mathbf{x}) &= B / A + p(\mathbf{x}) \\
   w(0,\mathbf{x}) &= 3.0 + p(\mathbf{x})

where the perturbation function is

.. math::

   p(\mathbf{x}) = \alpha e^{-((\mathbf{x}-\mu)^T \sigma^{-1}(\mathbf{x}-\mu)) /
   (2 \sqrt{|\sigma| 8 \pi^3}) }

with :math:`\alpha = 0.1`, :math:`\mu = 0.5\, \textbf{x}_\text{max}`, and
:math:`\sigma` is a diagonal matrix with entries
:math:`0.25\, \mathbf{x}_\text{max}`.

Spatial derivatives are discretized with first-order upwind finite differences
on a uniform spatial grid. The system can be evolved in time using explicit,
implicit, or IMEX methods from ARKODE, Adams or BDF methods from CVODE, or BDF
methods from IDA. When using an IMEX method, advection is treated explicitly and
reactions implicitly.

The nonlinear system(s) that arise in each time step may be solved using a
global Newton method with a matrix-free GMRES linear solver or an Anderson
accelerated fixed-point method. When using an IMEX method, a custom task-local
nonlinear solver that leverages the locality of the reaction systems may also be
used.


Options
^^^^^^^

Several command line options are available to change the problem parameters
as well as the integrator and solver options. A summary of the options are
listed in :numref:`Table.3D_advection_reaction_options`.

.. tabularcolumns:: |\Y{0.32}|\Y{0.53}|\Y{0.15}|

.. _Table.3D_advection_reaction_options:

.. Table:: 3D Advection-Reaction Benchmark Command Line Options

   +-------------------------------+------------------------------------+---------------+
   | Option                        | Description                        | Default       |
   +===============================+====================================+===============+
   | ``--help``                    | Print the command line options     | --            |
   |                               | and description                    |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--monitor``                 | Print solution information to      | Off           |
   |                               | the screen (slower)                |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--nout <int>``              | Number of output times             | 40            |
   +-------------------------------+------------------------------------+---------------+
   | ``--output-dir <dir>``        | Directory where all output files   | ``.``         |
   |                               | will be written                    |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--npts <int>``              | Number of mesh points in each      | 100           |
   |                               | direction                          |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--npxyz <int> <int> <int>`` | Number of MPI tasks in each        | 0 0 0         |
   |                               | direction (0 forces MPI to decide) |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--xmax <sunrealtype>``      | Maximum value of :math:`x`,        | 1.0           |
   |                               | :math:`y`, and :math:`z` in        |               |
   |                               | :math:`\textbf{x}_{\text{max}}`    |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--A <sunrealtype>``         | Constant concentration of species  | 1.0           |
   |                               | :math:`A`                          |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--B <sunrealtype>``         | Constant concentration of species  | 3.5           |
   |                               | :math:`B`                          |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--c <sunrealtype>``         | Advection speed :math:`c`          | 0.01          |
   +-------------------------------+------------------------------------+---------------+
   | ``--order <int>``             | Integration method order           | 3             |
   +-------------------------------+------------------------------------+---------------+
   | ``--method <method>``         | Integrator to use: ``ERK``,        | ``ARK-IMEX``  |
   |                               | ``ARK-DIRK``, ``ARK-IMEX``,        |               |
   |                               | ``CV-BDF``, ``CV-ADAMS``, ``IDA``  |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--fpaccel <int>``           | Number of fixed point acceleration | 3             |
   |                               | vectors                            |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--nls <method>``            | Nonlinear Solver Method:           | ``tl-newton`` |
   |                               | ``newton``, ``tl-newton``,         |               |
   |                               | ``fixedpoint``, ``none``           |               |
   +-------------------------------+------------------------------------+---------------+
   | ``--fused``                   | Enabled fused operations           | Off           |
   +-------------------------------+------------------------------------+---------------+
   | ``--tf <sunrealtype>``        | Final integration time :math:`t_f` | 10.0          |
   +-------------------------------+------------------------------------+---------------+
   | ``--rtol <sunrealtype>``      | Relative tolerance                 | 1.0e-6        |
   +-------------------------------+------------------------------------+---------------+
   | ``--atol <sunrealtype>``      | Absolute tolerance                 | 1.0e-9        |
   +-------------------------------+------------------------------------+---------------+


Building and Running
^^^^^^^^^^^^^^^^^^^^

To build the benchmark executables SUNDIALS must be configured with ARKODE,
CVODE, and IDA enabled and with MPI and RAJA support on. Additionally, either
CUDA or HIP support must be on to build executables utilizing NVIDIA or AMD
GPUs. See the installation guide for more details on configuring, building,
and installing SUNDIALS.

Based on the configuration the following executables will be built and installed
in the ``<install prefix>/bin/benchmarks/advection_reaction_3D`` directory:

* ``advection_reaction_3D`` -- MPI parallelism
* ``advection_reaction_3D_mpicuda`` -- MPI + CUDA parallelism
* ``advection_reaction_3D_mpihip`` -- MPI + HIP parallelism

On Summit, with the default environment

* Compiler: xl/16.1.1-5
* MPI: spectrum-mpi/10.3.1.2-20200121
* CUDA: cuda/10.1.243

an example ``jsrun`` command is

.. code-block:: none

   jsrun -n 2 -a 1 -c 1 -g 1 ./advection_reaction_3D_mpicuda

On Lassen, with the environment

* Compiler: gcc/8.3.1
* MPI: mvapich2/2021.05.28-cuda-11.1.1
* CUDA: cuda/11.1.1

an example ``jsrun`` command is

.. code-block:: none

   jsrun -n 2 -a 1 -c 1 -g 1 ./advection_reaction_3D_mpicuda
