.. ----------------------------------------------------------------
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

.. _ARKODE.Usage.PreconditionerModules:

Preconditioner modules
============================

The efficiency of Krylov iterative methods for the solution of linear
systems can be greatly enhanced through preconditioning.  For problems
in which the user cannot define a more effective, problem-specific
preconditioner, ARKODE provides two internal preconditioner modules:
a banded preconditioner for serial and threaded problems (ARKBANDPRE)
and a band-block-diagonal preconditioner for parallel problems (ARKBBDPRE).


.. _ARKODE.Usage.BandPre:

A serial banded preconditioner module
-------------------------------------------

This preconditioner provides a band matrix preconditioner for use with
iterative SUNLINSOL modules in a serial or threaded setting. It requires
that the problem be set up using either the
NVECTOR_SERIAL, NVECTOR_OPENMP or NVECTOR_PTHREADS module, due to data
access patterns.  It also currently requires that the problem involve
an identity mass matrix, i.e., :math:`M = I`.

This module uses difference quotients of the ODE right-hand
side function :math:`f^I` to generate a band matrix of bandwidth
``ml + mu + 1``, where the number of super-diagonals (``mu``, the
upper half-bandwidth) and sub-diagonals (``ml``, the lower
half-bandwidth) are specified by the user.  This band matrix is used
to to form a preconditioner the Krylov linear solver.  Although this
matrix is intended to approximate the Jacobian
:math:`J = \dfrac{\partial f^I}{\partial y}`, it may be a very crude
approximation, since the true Jacobian may not be banded, or its true
bandwidth may be larger than ``ml + mu + 1``.  However, as long as the
banded approximation generated for the preconditioner is sufficiently
accurate, it may speed convergence of the Krylov iteration.



ARKBANDPRE usage
"""""""""""""""""""""

In order to use the ARKBANDPRE module, the user need not define
any additional functions.  In addition to the header files required
for the integration of the ODE problem (see
:numref:`ARKODE.Usage.Headers`), to use the ARKBANDPRE module, the user's
program must include the header file ``arkode_bandpre.h`` which
declares the needed function prototypes.  The following is a summary
of the usage of this module.  Steps that are unchanged from the
skeleton program presented in :numref:`ARKODE.Usage.Skeleton` are
*italicized*.

#. *Initialize multi-threaded environment (if appropriate)*

#. *Create the SUNDIALS simulation context object.*

#. *Set problem dimensions*

#. *Set vector of initial values*

#. *Create ARKODE object*

#. *Specify integration tolerances*

#. Create iterative linear solver object

   When creating the iterative linear solver object, specify the type
   of preconditioning (``SUN_PREC_LEFT`` or ``SUN_PREC_RIGHT``) to use.

#. *Set linear solver optional inputs*

#. *Attach linear solver module*

#. Initialize the ARKBANDPRE preconditioner module

   Specify the upper and lower half-bandwidths (``mu`` and ``ml``,
   respectively) and call

   ``ier = ARKBandPrecInit(arkode_mem, N, mu, ml);``

   to allocate memory and initialize the internal preconditioner
   data.

#. *Create nonlinear solver object*

#. *Attach nonlinear solver module*

#. *Set nonlinear solver optional inputs*

#. *Set optional inputs*

   Note that the user should not call
   :c:func:`ARKodeSetPreconditioner()` as it will overwrite the
   preconditioner setup and solve functions.

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. Get optional outputs

   Additional optional outputs associated with ARKBANDPRE are
   available by way of the two routines described below,
   :c:func:`ARKBandPrecGetWorkSpace()` and
   :c:func:`ARKBandPrecGetNumRhsEvals()`.

#. *Deallocate memory for solution vector*

#. *Free solver memory*

#. *Free linear solver memory*

#. *Free nonlinear solver memory*





ARKBANDPRE user-callable functions
"""""""""""""""""""""""""""""""""""""

The ARKBANDPRE preconditioner module is initialized and attached
by calling the following function:



.. c:function:: int ARKBandPrecInit(void* arkode_mem, sunindextype N, sunindextype mu, sunindextype ml)

   Initializes the ARKBANDPRE preconditioner and
   allocates required (internal) memory for it.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param N: problem dimension (size of ODE system).
   :param mu: upper half-bandwidth of the Jacobian approximation.
   :param ml: lower half-bandwidth of the Jacobian approximation.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARKLS_MEM_FAIL: a memory allocation request failed.

   .. note::

      The banded approximate Jacobian will have nonzero elements
      only in locations :math:`(i,j)` with *ml* :math:`\le j-i \le` *mu*.



The following two optional output functions are available for use with
the ARKBANDPRE module:



.. c:function:: int ARKBandPrecGetWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the sizes of the ARKBANDPRE real and integer
   workspaces.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lenrwLS: the number of ``sunrealtype`` values in the
                   ARKBANDPRE workspace.
   :param leniwLS: the number of integer values in the  ARKBANDPRE workspace.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_PMEM_NULL: the preconditioner memory was ``NULL``.

   .. note::

      The workspace requirements reported by this routine
      correspond only to memory allocated within the ARKBANDPRE module
      (the banded matrix approximation, banded ``SUNLinearSolver``
      object, and temporary vectors).

      The workspaces referred to here exist in addition to those given by
      the corresponding function :c:func:`ARKodeGetLinWorkSpace()`.

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.



.. c:function:: int ARKBandPrecGetNumRhsEvals(void* arkode_mem, long int* nfevalsBP)

   Returns the number of calls made to the user-supplied
   right-hand side function :math:`f^I` for constructing the
   finite-difference banded Jacobian approximation used within the
   preconditioner setup function.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param nfevalsBP: number of calls to :math:`f^I`.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_PMEM_NULL: the preconditioner memory was ``NULL``.

   .. note::

      The counter *nfevalsBP* is distinct from the counter
      *nfevalsLS* returned by the corresponding function
      :c:func:`ARKodeGetNumLinRhsEvals()` and also from the number of
      evaluations returned by the time-stepping module (e.g., *nfi_evals*
      returned by :c:func:`ARKStepGetNumRhsEvals()`).  The total number of
      right-hand side function evaluations is the sum of all three of these
      counters.





.. _ARKODE.Usage.BBDPre:

A parallel band-block-diagonal preconditioner module
---------------------------------------------------------

A principal reason for using a parallel ODE solver (such as ARKODE)
lies in the solution of partial differential equations
(PDEs). Moreover, Krylov iterative methods are used on many such
problems due to the nature of the underlying linear system of
equations that needs to solved at each time step.  For many PDEs, the
linear algebraic system is large, sparse and structured.  However, if
a Krylov iterative method is to be effective in this setting, then a
nontrivial preconditioner is required.  Otherwise, the rate of
convergence of the Krylov iterative method is usually slow, and
degrades as the PDE mesh is refined.  Typically, an effective
preconditioner must be problem-specific.

However, we have developed one type of preconditioner that treats a
rather broad class of PDE-based problems.  It has been successfully
used with CVODE for several realistic, large-scale problems :cite:p:`HiTa:98`,
and is included in a software module within the ARKODE package.  This
preconditioning module works with the parallel vector module
NVECTOR_PARALLEL and is usable with any of the Krylov iterative linear
solvers through the ARKLS interface. It generates a preconditioner
that is a block-diagonal matrix with each block being a band
matrix. The blocks need not have the same number of super- and
sub-diagonals and these numbers may vary from block to block. This
Band-Block-Diagonal Preconditioner module is called ARKBBDPRE.

One way to envision these preconditioners is to think of the
computational PDE domain as being subdivided into :math:`Q`
non-overlapping subdomains, where each subdomain is assigned to one of
the :math:`Q` MPI tasks used to solve the ODE system.  The basic idea
is to isolate the preconditioning so that it is local to each process,
and also to use a (possibly cheaper) approximate right-hand side
function for construction of this preconditioning matrix.  This
requires the definition of a new function :math:`g(t,y) \approx
f^I(t,y)` that will be used to construct the BBD preconditioner
matrix.  At present, we assume that the ODE be written in explicit
form as

.. math::
   \dot{y} = f^E(t,y) + f^I(t,y),

where :math:`f^I` corresponds to the ODE components to be treated
implicitly, i.e. this preconditioning module does not support problems
with non-identity mass matrices.  The user may set :math:`g = f^I`, if
no less expensive approximation is desired.

Corresponding to the domain decomposition, there is a decomposition of
the solution vector :math:`y` into :math:`Q` disjoint blocks
:math:`y_q`, and a decomposition of :math:`g` into blocks
:math:`g_q`. The block :math:`g_q` depends both on :math:`y_p` and on
components of blocks :math:`y_{q'}` associated with neighboring
subdomains (so-called ghost-cell data).  If we let :math:`\bar{y}_q`
denote :math:`y_q` augmented with those other components on which
:math:`g_q` depends, then we have

.. math::
   g(t,y) = \left[ g_1(t,\bar{y}_1), g_2(t,\bar{y}_2), \ldots , g_Q(t,\bar{y}_Q) \right]^T,

and each of the blocks :math:`g_q(t,\bar{y}_q)` is decoupled from one another.

The preconditioner associated with this decomposition has the form

.. math::
   P = \begin{bmatrix} P_1 & & & \\ & P_2 & & \\ & & \ddots &\\ & & & P_Q \end{bmatrix}

where

.. math::
   P_q \approx I - \gamma J_q

and where :math:`J_q` is a difference quotient approximation to
:math:`\dfrac{\partial g_q}{\partial \bar{y}_q}`.  This matrix is taken
to be banded, with upper and lower half-bandwidths *mudq* and
*mldq* defined as the number of non-zero diagonals above and below
the main diagonal, respectively.  The difference quotient
approximation is computed using *mudq* + *mldq* + 2 evaluations of
:math:`g_m`, but only a matrix of bandwidth *mukeep* + *mlkeep* + 1 is
retained. Neither pair of parameters need be the true half-bandwidths
of the Jacobian of the local block of :math:`g`, if smaller values
provide a more efficient preconditioner. The solution of the complete
linear system

.. math::
   Px = b

reduces to solving each of the distinct equations

.. math::
   P_q x_q = b_q, \quad q=1,\ldots,Q,

and this is done by banded LU factorization of :math:`P_q` followed by
a banded backsolve.

Similar block-diagonal preconditioners could be considered with
different treatments of the blocks :math:`P_q`.  For example,
incomplete LU factorization or an iterative method could be used
instead of banded LU factorization.



ARKBBDPRE user-supplied functions
""""""""""""""""""""""""""""""""""

The ARKBBDPRE module calls two user-provided functions to construct
:math:`P`: a required function *gloc* (of type :c:func:`ARKLocalFn()`)
which approximates the right-hand side function :math:`g(t,y) \approx
f^I(t,y)` and which is computed locally, and an optional function
*cfn* (of type :c:func:`ARKCommFn()`) which performs all inter-process
communication necessary to evaluate the approximate right-hand side
:math:`g`. These are in addition to the user-supplied right-hand side
function :math:`f^I`. Both functions take as input the same pointer
*user_data* that is passed by the user to
:c:func:`ARKodeSetUserData()` and that was passed to the user's
function :math:`f^I`. The user is responsible for providing space
(presumably within *user_data*) for components of :math:`y` that are
communicated between processes by *cfn*, and that are then used by
*gloc*, which should not do any communication.



.. c:type:: int (*ARKLocalFn)(sunindextype Nlocal, sunrealtype t, N_Vector y, N_Vector glocal, void* user_data)

   This *gloc* function computes :math:`g(t,y)`.  It
   fills the vector *glocal* as a function of *t* and *y*.

   :param Nlocal: the local vector length.
   :param t: the value of the independent variable.
   :param y: the value of the dependent variable vector on this process.
   :param glocal: the output vector of :math:`g(t,y)` on this process.
   :param user_data: a pointer to user data, the same as the
                     *user_data* parameter passed to :c:func:`ARKodeSetUserData()`.

   :return: An *ARKLocalFn* should return 0 if successful, a positive value if
            a recoverable error occurred (in which case ARKODE will attempt to
            correct), or a negative value if it failed unrecoverably (in which
            case the integration is halted and :c:func:`ARKodeEvolve()` will return
            *ARK_LSETUP_FAIL*).

   .. note::

      This function should assume that all inter-process
      communication of data needed to calculate *glocal* has already been
      done, and that this data is accessible within user data.

      The case where :math:`g` is mathematically identical to :math:`f^I`
      is allowed.



.. c:type:: int (*ARKCommFn)(sunindextype Nlocal, sunrealtype t, N_Vector y, void* user_data)

   This *cfn* function performs all inter-process
   communication necessary for the execution of the *gloc* function
   above, using the input vector *y*.

   :param Nlocal: the local vector length.
   :param t: the value of the independent variable.
   :param y: the value of the dependent variable vector on this process.
   :param user_data: a pointer to user data, the same as the
                     *user_data* parameter passed to :c:func:`ARKodeSetUserData()`.

   :return: An *ARKCommFn* should return 0 if successful, a positive value if a
            recoverable error occurred (in which case ARKODE will attempt to
            correct), or a negative value if it failed unrecoverably (in which
            case the integration is halted and :c:func:`ARKodeEvolve()` will return
            *ARK_LSETUP_FAIL*).

   .. note::

      The *cfn* function is expected to save communicated data in
      space defined within the data structure *user_data*.

      Each call to the *cfn* function is preceded by a call to the
      right-hand side function :math:`f^I` with the same :math:`(t,y)`
      arguments. Thus, *cfn* can omit any communication done by
      :math:`f^I` if relevant to the evaluation of *glocal*. If all
      necessary communication was done in :math:`f^I`, then *cfn* =
      ``NULL`` can be passed in the call to :c:func:`ARKBBDPrecInit()`
      (see below).




ARKBBDPRE usage
"""""""""""""""""""""

In addition to the header files required for the integration of the
ODE problem (see :numref:`ARKODE.Usage.Headers`), to use the
ARKBBDPRE module, the user's program must include the header file
``arkode_bbdpre.h`` which declares the needed function prototypes.

The following is a summary of the proper usage of this module. Steps
that are unchanged from the skeleton program presented in
:numref:`ARKODE.Usage.Skeleton` are *italicized*.

#. *Initialize MPI*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions*

#. *Set vector of initial values*

#. *Create ARKODE object*

#. *Specify integration tolerances*

#. Create iterative linear solver object

   When creating the iterative linear solver object, specify the type
   of preconditioning (``SUN_PREC_LEFT`` or ``SUN_PREC_RIGHT``) to use.

#. *Set linear solver optional inputs*

#. *Attach linear solver module*

#. Initialize the ARKBBDPRE preconditioner module

   Specify the upper and lower half-bandwidths for computation
   ``mudq`` and ``mldq``, the upper and lower half-bandwidths for
   storage ``mukeep`` and ``mlkeep``, and call

   ``ier = ARKBBDPrecInit(arkode_mem, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn);``

   to allocate memory and initialize the internal preconditioner
   data. The last two arguments of :c:func:`ARKBBDPrecInit()` are the
   two user-supplied functions of type :c:func:`ARKLocalFn()` and
   :c:func:`ARKCommFn()` described above, respectively.

#. *Create nonlinear solver object*

#. *Attach nonlinear solver module*

#. *Set nonlinear solver optional inputs*

#. *Set optional inputs*

   Note that the user should not call
   :c:func:`ARKodeSetPreconditioner()` as it will overwrite the
   preconditioner setup and solve functions.

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. *Get optional outputs*

   Additional optional outputs associated with ARKBBDPRE are
   available through the routines
   :c:func:`ARKBBDPrecGetWorkSpace()` and
   :c:func:`ARKBBDPrecGetNumGfnEvals()`.

#. *Deallocate memory for solution vector*

#. *Free solver memory*

#. *Free linear solver memory*

#. *Free nonlinear solver memory*

#. *Finalize MPI*





ARKBBDPRE user-callable functions
""""""""""""""""""""""""""""""""""""

The ARKBBDPRE preconditioner module is initialized (or re-initialized)
and attached to the integrator by calling the following functions:

.. c:function:: int ARKBBDPrecInit(void* arkode_mem, sunindextype Nlocal, sunindextype mudq, sunindextype mldq, sunindextype mukeep, sunindextype mlkeep, sunrealtype dqrely, ARKLocalFn gloc, ARKCommFn cfn)

   Initializes and allocates (internal) memory for the
   ARKBBDPRE preconditioner.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param Nlocal: local vector length.
   :param mudq: upper half-bandwidth to be used in the difference
                quotient Jacobian approximation.
   :param mldq: lower half-bandwidth to be used in the difference
                quotient Jacobian approximation.
   :param mukeep: upper half-bandwidth of the retained banded
                  approximate Jacobian block.
   :param mlkeep: lower half-bandwidth of the retained banded
                  approximate Jacobian block.
   :param dqrely: the relative increment in components of *y* used in
                  the difference quotient approximations.  The default is *dqrely*
                  = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
                  passing *dqrely* = 0.0.
   :param gloc: the name of the C function (of type :c:func:`ARKLocalFn()`)
                which computes the approximation :math:`g(t,y) \approx f^I(t,y)`.
   :param cfn: the name of the C function (of type :c:func:`ARKCommFn()`) which
               performs all inter-process communication required for the
               computation of :math:`g(t,y)`.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_ILL_INPUT: an input had an illegal value.
   :retval ARKLS_MEM_FAIL: a memory allocation request failed.

   .. note::

      If one of the half-bandwidths *mudq* or *mldq* to be used
      in the difference quotient calculation of the approximate Jacobian is
      negative or exceeds the value *Nlocal*-1, it is replaced by 0 or
      *Nlocal*-1 accordingly.

      The half-bandwidths *mudq* and *mldq* need not be the true
      half-bandwidths of the Jacobian of the local block of :math:`g`
      when smaller values may provide a greater efficiency.

      Also, the half-bandwidths *mukeep* and *mlkeep* of the retained
      banded approximate Jacobian block may be even smaller than
      *mudq* and *mldq*, to reduce storage and computational costs
      further.

      For all four half-bandwidths, the values need not be the same on
      every processor.



The ARKBBDPRE module also provides a re-initialization function to
allow solving a sequence of problems of the same size, with the same
linear solver choice, provided there is no change in *Nlocal*,
*mukeep*, or *mlkeep*. After solving one problem, and after
calling ``*StepReInit`` to re-initialize ARKODE for a
subsequent problem, a call to :c:func:`ARKBBDPrecReInit()` can be made
to change any of the following: the half-bandwidths *mudq* and
*mldq* used in the difference-quotient Jacobian approximations, the
relative increment *dqrely*, or one of the user-supplied functions
*gloc* and *cfn*. If there is a change in any of the linear solver
inputs, an additional call to the "Set" routines provided by the
SUNLINSOL module, and/or one or more of the corresponding
``ARKodeSet***`` functions, must also be made (in the proper order).


.. c:function:: int ARKBBDPrecReInit(void* arkode_mem, sunindextype mudq, sunindextype mldq, sunrealtype dqrely)

   Re-initializes the ARKBBDPRE preconditioner module.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param mudq: upper half-bandwidth to be used in the difference
                quotient Jacobian approximation.
   :param mldq: lower half-bandwidth to be used in the difference
                quotient Jacobian approximation.
   :param dqrely: the relative increment in components of *y* used in
                  the difference quotient approximations.  The default is *dqrely*
                  = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
                  passing *dqrely* = 0.0.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_PMEM_NULL: the preconditioner memory was ``NULL``.

   .. note::

      If one of the half-bandwidths *mudq* or *mldq* is
      negative or exceeds the value *Nlocal*-1, it is replaced by 0 or
      *Nlocal*-1 accordingly.


The following two optional output functions are available for use with
the ARKBBDPRE module:


.. c:function:: int ARKBBDPrecGetWorkSpace(void* arkode_mem, long int* lenrwBBDP, long int* leniwBBDP)

   Returns the processor-local ARKBBDPRE real and
   integer workspace sizes.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param lenrwBBDP: the number of ``sunrealtype`` values in the
                     ARKBBDPRE workspace.
   :param leniwBBDP: the number of integer values in the  ARKBBDPRE workspace.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_PMEM_NULL: the preconditioner memory was ``NULL``.

   .. note::

      The workspace requirements reported by this routine
      correspond only to memory allocated within the ARKBBDPRE module
      (the banded matrix approximation, banded ``SUNLinearSolver``
      object, temporary vectors). These values are local to each process.

      The workspaces referred to here exist in addition to those given by
      the corresponding function :c:func:`ARKodeGetLinWorkSpace()`.

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.



.. c:function:: int ARKBBDPrecGetNumGfnEvals(void* arkode_mem, long int* ngevalsBBDP)

   Returns the number of calls made to the user-supplied
   *gloc* function (of type :c:func:`ARKLocalFn()`) due to the finite
   difference approximation of the Jacobian blocks used within the
   preconditioner setup function.

   :param arkode_mem: pointer to the ARKODE memory block.
   :param ngevalsBBDP: the number of calls made to the user-supplied
                       *gloc* function.

   :retval ARKLS_SUCCESS: the function exited successfully.
   :retval ARKLS_MEM_NULL: ``arkode_mem`` was ``NULL``.
   :retval ARKLS_LMEM_NULL: the linear solver memory was ``NULL``.
   :retval ARKLS_PMEM_NULL: the preconditioner memory was ``NULL``.


In addition to the *ngevalsBBDP* *gloc* evaluations, the costs
associated with ARKBBDPRE also include *nlinsetups* LU
factorizations, *nlinsetups* calls to *cfn*, *npsolves* banded
backsolve calls, and *nfevalsLS* right-hand side function
evaluations, where *nlinsetups* is an optional ARKODE output and
*npsolves* and *nfevalsLS* are linear solver optional outputs (see
the table :numref:`ARKODE.Usage.ARKLsOutputs`).
