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

:tocdepth: 3

.. _parhyp_c:

====================================
Parallel Hypre example problems
====================================




.. _ark_diurnal_kry_ph:

ark_diurnal_kry_ph
============================================

This problem is mathematically identical to the parallel C example
problem :ref:`ark_diurnal_kry_p`.  As before, this test problem models
a two-species diurnal kinetics advection-diffusion PDE system in two
spatial dimensions,

.. math::

   \frac{\partial c_i}{\partial t} =
     K_h \frac{\partial^2 c_i}{\partial x^2} +
     V \frac{\partial     c_i}{\partial x} +
     \frac{\partial}{\partial y}\left( K_v(y)
     \frac{\partial c_i}{\partial y}\right) +
     R_i(c_1,c_2,t),\quad i=1,2

where

.. math::

   R_1(c_1,c_2,t) &= -q_1*c_1*c_3 - q_2*c_1*c_2 + 2*q_3(t)*c_3 + q_4(t)*c_2, \\
   R_2(c_1,c_2,t) &=  q_1*c_1*c_3 - q_2*c_1*c_2 - q_4(t)*c_2, \\
   K_v(y) &= K_{v0} e^{y/5}.

Here :math:`K_h`, :math:`V`, :math:`K_{v0}`, :math:`q_1`, :math:`q_2`,
and :math:`c_3` are constants, and :math:`q_3(t)` and :math:`q_4(t)`
vary diurnally.  The problem is posed on the square spatial domain
:math:`(x,y) \in [0,20]\times[30,50]`, with homogeneous Neumann
boundary conditions, and for time interval :math:`t\in [0,86400]` sec
(1 day).

We enforce the initial conditions

.. math::

   c^1(x,y) &=  10^6 \chi(x)\eta(y) \\
   c^2(x,y) &=  10^{12} \chi(x)\eta(y) \\
   \chi(x) &= 1 - \sqrt{\frac{x - 10}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}} \\
   \eta(y) &= 1 - \sqrt{\frac{y - 40}{10}} + \frac12 \sqrt[4]{\frac{x - 10}{10}}.




Numerical method
----------------

The numerical method is identical to the previous implementation,
except that we now use the HYPRE parallel vector module,
NVECTOR_PARHYP.  The output of these two examples is identical.
Below, we discuss only the main differences between the two
implementations; familiarity with the HYPRE library is helpful.

We use the HYPRE IJ vector interface to allocate the template vector
and create the parallel partitioning:

.. code-block:: c

   HYPRE_IJVectorCreate(comm, my_pe*local_N, (my_pe + 1)*local_N - 1, &Uij);
   HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(Uij);

The *initialize* call means that vector elements are ready to be set using
the IJ interface. We choose the initial condition vector :math:`x_0 =
x(t_0)` as the template vector, and we set its values in the
``SetInitialProfiles(...)`` function. We complete assembly of the
HYPRE vector with the calls:

.. code-block:: c

   HYPRE_IJVectorAssemble(Uij);
   HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

The *assemble* call is collective and makes the HYPRE vector ready to
use.  The sets the handle ``Upar`` to the actual HYPRE vector.  The
handle is then passed to the ``N_VMake`` function, which creates the
template ``N_Vector``, ``u``, as a wrapper around the HYPRE vector.
All of the other vectors used in the computation are created by
cloning this template vector.

Furthermore, since the template vector does not own the underlying
HYPRE vector (it was created using the ``HYPRE_IJVectorCreate`` call
above), so it is the user's responsibility to destroy it by calling
``HYPRE_IJVectorDestroy(Uij)`` after the template vector ``u`` has
been destroyed.  This function will destroy both the HYPRE vector
and its IJ interface.

To access individual elements of the solution and derivative vectors
``u`` and ``udot`` in the IVP right-hand side function, ``f``, the
user needs to first extract the HYPRE vector by calling
``N_VGetVector_ParHyp``, and then use HYPRE-specific methods to access
the data from that point on.

.. note::

   Currently, interfaces to HYPRE solvers and preconditioners are not
   available directly through the SUNDIALS interfaces, however these
   could be utilized directly as preconditioners for SUNDIALS' Krylov
   solvers.  Direct interfaces to the HYPRE solvers will be provided
   in subsequent SUNDIALS releases.  The current HYPRE vector
   interface is included in this release mainly for testing purposes
   and as a preview of functionality to come.
