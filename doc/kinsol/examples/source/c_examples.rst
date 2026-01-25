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

.. _KINSOL.Examples.C:

C Example Problems
==================

A serial dense example: kinFerTron_dns
--------------------------------------

As an initial illustration of the use of the KINSOL package for the solution of
nonlinear systems, we give a sample program called ``kinFerTron_dns.c``. It uses
the :ref:`dense linear solver <SUNLinSol_Dense>` and the :ref:`serial vector
<NVectors.NVSerial>` for the solution of the Ferraris-Tronconi test problem
:cite:p:`FlPa:99`.

This problem involves a blend of trigonometric and exponential terms,

.. math::

   0 &= 0.5 \sin(x_1 x_2) - 0.25 x_2/\pi - 0.5 x_1, \\
   0 &= (1-0.25/\pi) ( e^{2 x_1} - e ) + e x_2 / \pi - 2 e x_1, \\

and is subject to the following constraints

.. math::

   x_{1\min} &= 0.25 \le x_1 \le 1 = x_{1\max} \\
   x_{2\min} &= 1.5 \le x_2 \le 2\pi = x_{2\max}

The bounds constraints on :math:`x_1` and :math:`x_2` are treated by introducing
four additional variables and using KINSOL's optional constraints feature to
enforce non-positivity and non-negativity:

.. math::

   l_1 &= x_1 - x_{1\min} \ge 0\\
   L_1 &= x_1 - x_{1\max} \le 0\\
   l_2 &= x_2 - x_{2\min} \ge 0\\
   L_2 &= x_2 - x_{2\max} \le 0

The Ferraris-Tronconi problem has two known solutions. We solve it with KINSOL
using two sets of initial guesses for :math:`x_1` and :math:`x_2` (first their
lower bounds and secondly the middle of their feasible regions), both with an
exact and modified Newton method, with and without line search.

Following the initial comment block, this program has a number of ``#include``
lines, which allow access to useful items in KINSOL header files:

* ``kinsol.h`` provides prototypes for the KINSOL functions to be called and
  also a number of constants that are to be used in setting input arguments and
  testing the return value of :c:func:`KINSol` .

* ``nvector_serial.h`` provides prototypes for the serial vector implementation.

* ``sunmatrix_dense.h`` provides the prototypes for the dense matrix
  implementation.

* ``sunlinsol_dense.h`` provides the prototypes for the dense linear solver
  implementation.

* ``sundials_types.h`` file provides the definition of :c:type:`sunrealtype`.
  For now, it suffices to read :c:type:`sunrealtype` as ``double``.

Next, the program defines some problem-specific constants, which are isolated to
this early location to make it easy to change them as needed.

The program prologue ends with prototypes of the user-supplied system function
``func`` and several helper functions.

The ``main`` program begins with creating the SUNDIALS context object
(``sunctx``) which must be passed to all other SUNDIALS constructors. Then we
allocated the user data structure, ``data``, which contains two arrays with
lower and upper bounds for :math:`x_1` and :math:`x_2`. Next, we create five
serial vectors for the two different initial guesses (``u1`` and ``u2``), the
solution vector (``u``), the scaling factors (``s``), and the constraint
specifications (``c``).

The initial guess vectors ``u1`` and ``u2`` are filled by the functions
``SetInitialGuess1`` and ``SetInitialGuess2`` and the constraint vector ``c`` is
initialized to ``[0,0,1,-1,1,-1]`` indicating that there are no additional
constraints on the first two components of ``u`` (i.e., :math:`x_1` and
:math:`x_2`) and that the 3rd and 5th components should be non-negative, while
the 4th and 6th should be non-positive.

The calls to :c:func:`N_VNew_Serial`, and also later calls to various ``KIN***``
functions, make use of a function, ``check_flag``, which examines the
return value and prints a message if there was a failure.  The ``check_flag``
function was written to be used for any serial SUNDIALS application.

The call to :c:func:`KINCreate` creates the KINSOL solver memory block. Its
return value is a pointer to that memory block for this problem. In the case of
failure, the return value is ``NULL``. This pointer must be passed in the
remaining calls to KINSOL functions.

The next four calls to KINSOL optional input functions specify the pointer to
the user data structure (to be passed to all subsequent calls to ``func``), the
vector of additional constraints, and the function and scaled step tolerances,
``fnormtol`` and ``scsteptol``, respectively.

The solver is initialized with :c:func:`KINInit` call which specifies the system
function ``func`` and provides the vector ``u`` which will be used internally as
a template for cloning additional necessary vectors of the same type as ``u``.

The call to :c:func:`SUNDenseMatrix` to creates the Jacobian matrix to use with
the dense linear solver which is created by :c:func:`SUNLinSol_Dense`. Finally,
both of these objects are attached to KINSOL by calling
:c:func:`KINSetLinearSolver`.

The main program proceeds by solving the nonlinear system eight times, using
each of the two initial guesses, ``u1`` and ``u2`` (which are first copied into
the vector ``u`` using :c:func:`N_VScale`), with and without globalization
through line search (specified by setting ``glstr`` to ``KIN_LINESEARCH`` and
``KIN_NONE``, respectively), and applying either an exact or a modified Newton
method. The switch from exact to modified Newton is done by changing the number
of nonlinear iterations after which a Jacobian evaluation is enforced, a value
``mset=1`` thus resulting in re-evaluating the Jacobian at every single
iteration of the nonlinear solver (exact Newton method). Note that passing
``mset=0`` indicates using the default KINSOL value of 10.

The actual problem solution is carried out in the function ``SolveIt`` which
calls the main solver function, :c:func:`KINSol`, after first setting the
optional input ``mset``. After a successful return from :c:func:`KINSol`, the
solution :math:`[x_1, x_2]` and some solver statistics are printed.

The function ``func`` is a straightforward expression of the extended nonlinear
system. It uses the :c:func:`N_VGetArrayPointer` function to extract the data
arrays of the vectors ``u`` and ``f`` and sets the components of ``fdata`` using
the current values for the components of ``udata``. See :c:type:`KINSysFn` for a
detailed specification of ``func``.

The output generated by ``kinFerTron_dns`` is shown below.

.. literalinclude:: ../../../../examples/kinsol/serial/kinFerTron_dns.out
   :language: none

A serial Krylov example: kinFoodWeb_kry
---------------------------------------

We give here an example that illustrates the use of KINSOL with the GMRES Krylov
method as the linear system solver.

This program solves a nonlinear system that arises from a discretized system of
partial differential equations. The PDE system is a six-species food web
population model, with predator-prey interaction and diffusion on the unit
square in two dimensions. Given the dependent variable vector of species
concentrations :math:`c = [c_1, c_2,..., c_{n_s}]^T`, where :math:`n_s = 2 n_p`
is the number of species and :math:`n_p` is the number of predators and of prey,
then the PDEs can be written as

.. math::

   d_i \cdot \left( \frac{\partial^2 c_i}{\partial x^2} +
     \frac{\partial^2 c_i}{\partial y^2} \right) + f_i(x,y,c) = 0
   \quad (i=1,...,n_s)

where the subscripts :math:`i` are used to distinguish the species, and where

.. math::

   f_i(x,y,c) = c_i \cdot \left(b_i + \sum_{j=1}^{n_s} a_{i,j} \cdot c_j \right)

The problem coefficients are given by

.. math::

   a_{ij} =
   \begin{cases}
     -1                 & i=j \\
     -0.5 \cdot 10^{-6} & i \leq n_p , ~ j > n_p  \\
     10^4               & i > n_p , ~ j \leq n_p  \\
     0                  & \text{all other}
   \end{cases}

.. math::

   b_i = b_i(x,y) =
   \begin{cases}
     1 + \alpha xy   & i \leq n_p  \\
     -1 - \alpha xy   & i > n_p
   \end{cases}

and

.. math::

   d_i =
   \begin{cases}
     1 & i \leq n_p  \\
     0.5 & i > n_p
   \end{cases}

The spatial domain is the unit square :math:`(x,y) \in [0,1] \times [0,1]`.

Homogeneous Neumann boundary conditions are imposed and the initial guess is
constant in both :math:`x` and :math:`y`. For this example, the equations are
discretized spatially with standard central finite differences on a :math:`8
\times 8` mesh with :math:`n_s = 6`, giving a system of size 384.

Among the initial ``#include`` lines in this case is ``sunlinsol_spgmr.h`` which
contains constants and function prototypes associated with the :ref:`GMRES
linear solver <SUNLinSol.SPGMR>`.

The ``main`` program calls :c:func:`KINCreate` and then calls :c:func:`KINInit`
with the name of the user-supplied system function ``func`` and solution vector
as arguments.  The ``main`` program then calls a number of ``KINSet*`` routines
to notify KINSOL of the user data pointer, the positivity constraints on the
solution, and convergence tolerances on the system function and step size. It
calls :c:func:`SUNLinSol_SPGMR` to create the linear solver, supplying the
``maxl`` value of 15 as the maximum Krylov subspace dimension. It then calls
:c:func:`KINSetLinearSolver` to attach this solver module to KINSOL.  Next, a
maximum value of ``maxlrst = 2`` restarts is imposed through a call to
:c:func:`SUNLinSol_SPGMRSetMaxRestarts`. Finally, the user-supplied
preconditioner setup and solve functions, ``PrecSetupBD`` and ``PrecSolveBD``,
are specified through a call to :c:func:`KINSetPreconditioner`.  The ``data``
pointer passed to :c:func:`KINSetUserData` is passed to ``PrecSetupBD`` and
``PrecSolveBD`` whenever these are called.

Next, :c:func:`KINSol` is called to solve the system, the return value is tested
for error conditions, and the approximate solution vector is printed via a call
to ``PrintOutput``. After that, ``PrintFinalStats`` is called to get and print
final statistics, and memory is freed by calls to :c:func:`N_VDestroy`,
``FreeUserData``, :c:func:`KINFree`, and :c:func:`SUNContext_Free`. The
statistics printed are the total numbers of nonlinear iterations (``nni``),
``func`` evaluations (``nfe``, excluding those for :math:`Jv` product
evaluations), ``func`` evaluations for :math:`Jv` evaluations (``nfeSG``),
linear (Krylov) iterations (``nli``), preconditioner evaluations (``npe``), and
preconditioner solves (``nps``).  See the :ref:`KINSOL.Usage.CC.optional_output`
section for more information.

Mathematically, the dependent variable has three dimensions: species number,
:math:`x` mesh point, and :math:`y` mesh point. The macro ``IJ_Vptr`` isolates
the translation from three dimensions to the one-dimensional contiguous array of
components under the serial vector. This macro allows for clearer code and makes
it easy to change the underlying layout of the three-dimensional data.

The preconditioner used here is the block-diagonal part of the true Newton
matrix and is based only on the partial derivatives of the interaction terms
:math:`f` in the above equation and hence its diagonal blocks are :math:`n_s
\times n_s` matrices (:math:`n_s = 6`). It is generated and factored in the
``PrecSetupBD`` routine and backsolved in the ``PrecSolveBD`` routine. See
:c:type:`KINLsPrecSetupFn` and :c:type:`KINLsPrecSolveFn` for detailed
descriptions of these preconditioner functions.

The program ``kinFoodWeb_kry.c`` uses the "small" dense functions for all
operations on the :math:`6 \times 6` preconditioner blocks. Thus it includes
``sundials_dense.h``, and calls the small dense matrix functions
``SUNDlsMat_newDenseMat``, ``SUNDlsMat_denseGETRF``, and
``SUNDlsMat_denseGETRS``. The small dense functions are generally available for
KINSOL user programs (for more information, see the comments in the header file
``sundials_dense.h``).

In addition to the functions called by KINSOL, ``kinFoodWeb_kry.c`` includes
definitions of several functions. These are: ``AllocUserData`` to allocate
space for the preconditioner and the pivot arrays; ``InitUserData`` to load
problem constants in the ``data`` block; ``FreeUserData`` to free that block;
``SetInitialProfiles`` to load the initial values in ``cc``; ``PrintHeader`` to
print the heading for the output; ``PrintOutput`` to retrieve and print selected
solution values; ``PrintFinalStats`` to print statistics; and ``check_flag`` to
check return values for error conditions.

The output generated by ``kinFoodWeb_kry`` is shown below.  Note that the
solution involved 9 Newton iterations, with an average of about 37 Krylov
iterations per Newton iteration.

.. literalinclude:: ../../../../examples/kinsol/serial/kinFoodWeb_kry.out
   :language: none

A parallel example: kinFoodWeb_kry_bbd_p
----------------------------------------

In this example, ``kinFoodWeb_kry_bbd_p``, we solve the same problem as with
``kinFoodWeb_kry`` above, but in parallel, and instead of supplying the
preconditioner we use the :ref:`KINBBDPRE <KINSOL.Usage.CC.kin_bbdpre>`
preconditioner.

In this case, we think of the parallel MPI processes as being laid out in a
rectangle, and each process being assigned a subgrid of size ``MXSUB``
:math:`\times` ``MYSUB`` of the x-y grid. If there are ``NPEX`` processes in the
x direction and ``NPEY`` processes in the y direction, then the overall grid
size is ``MX`` :math:`\times` ``MY`` with ``MX = NPEX * MXSUB`` and ``MY =
NPEY * MYSUB``, and the size of the nonlinear system is ``NUM_SPECIES * MX *
MY``.

The evaluation of the nonlinear system function is performed in ``func``.  In
this parallel setting, the processes first communicate the subgrid boundary data
and then compute the local components of the nonlinear system function. The MPI
communication is isolated in the function ``ccomm`` (which in turn calls
``BRecvPost``, ``BSend``, and ``BRecvWait``) and the subgrid boundary data
received from neighboring processes is loaded into the work array ``cext``. The
computation of the nonlinear system function is done in ``func_local`` which
starts by copying the local segment of the ``cc`` vector into ``cext``, and then
by imposing the boundary conditions by copying the first interior mesh line from
``cc`` into ``cext``. After this, the nonlinear system function is evaluated by
using central finite-difference approximations using the data in ``cext``
exclusively.

:ref:`KINBBDPRE <KINSOL.Usage.CC.kin_bbdpre>` uses a band-block-diagonal
preconditioner, generated by difference quotients. The upper and lower
half-bandwidths of the Jacobian block generated on each process are both equal
to :math:`2 n_s - 1`, and that is the value passed as ``mudq`` and ``mldq`` in
the call to :c:func:`KINBBDPrecInit`.  These values are much less than the true
half-bandwidths of the Jacobian blocks, which are :math:`n_s \times` ``MXSUB``.
However, an even narrower band matrix is retained as the preconditioner, with
half-bandwidths equal to :math:`n_s`, and this is the value passed to
:c:func:`KINBBDPrecInit` for ``mukeep`` and ``mlkeep``.

The function ``func_local`` is also passed as the ``gloc`` argument to
:c:func:`KINBBDPrecInit`. Since all communication needed for the evaluation of
the local approximation of :math:`f` used in building the band-block-diagonal
preconditioner is already done for the evaluation of :math:`f` in ``func``, a
``NULL`` pointer is passed as the ``gcomm`` argument to
:c:func:`KINBBDPrecInit`.

The ``main`` program resembles closely that of the ``kinFoodWeb_kry`` example,
with particularization arising from the use of the :ref:`MPI parallel vector
<NVectors.NVParallel>`. It begins by initializing MPI and obtaining the total
number of processes and the rank of the local process. The local length of the
solution vector is then computed as ``NUM_SPECIES * MXSUB * MYSUB``. Distributed
vectors are created by calling the constructor :c:func:`N_VNew_Parallel` with
the MPI communicator and the local and global problem sizes as arguments. All
output is performed only from the process with ID equal to 0. Finally, after
:c:func:`KINSol` is called and the results are printed, all memory is
deallocated, and the MPI environment is terminated by calling ``MPI_Finalize``.

The output generated by ``kinFoodWeb_kry_bbd_p`` is shown below. Note that 9
Newton iterations were required, with an average of about 52 Krylov iterations
per Newton iteration.

.. literalinclude:: ../../../../examples/kinsol/parallel/kinFoodWeb_kry_bbd_p.out
   :language: none
