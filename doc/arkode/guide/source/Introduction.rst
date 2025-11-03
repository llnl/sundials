.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ UMBC
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

.. _ARKODE.Introduction:

************
Introduction
************

The ARKODE infrastructure provides adaptive-step time integration
modules for stiff, nonstiff and mixed stiff/nonstiff systems of
ordinary differential equations (ODEs).  ARKODE itself is structured
to support a wide range of one-step (but multi-stage) methods,
allowing for rapid development of parallel implementations of
state-of-the-art time integration methods.  At present, ARKODE is
packaged with four time-stepping modules, *ARKStep*, *ERKStep*, *SPRKStep*,
and *MRIStep*.


*ARKStep* supports ODE systems posed in split, linearly-implicit form,

.. math::
   M(t)\, \dot{y} = f^E(t,y) + f^I(t,y),  \qquad y(t_0) = y_0,
   :label: ARKODE_ODE_split_linearly_implicit

where :math:`t` is the independent variable, :math:`y` is the set of
dependent variables (in :math:`\mathbb{R}^N`), :math:`M` is a
user-specified, nonsingular operator from :math:`\mathbb{R}^N` to
:math:`\mathbb{R}^N`, and the right-hand side function is partitioned
into up to two components:

- :math:`f^E(t,y)` contains the "nonstiff" time scale components to be
  integrated explicitly, and
- :math:`f^I(t,y)`  contains the "stiff" time scale components to be
  integrated implicitly.

Either of these operators may be disabled, allowing for fully
explicit, fully implicit, or combination implicit-explicit (ImEx) time
integration.

The algorithms used in ARKStep are adaptive- and fixed-step additive
Runge--Kutta methods. Such methods are defined through combining two
complementary Runge--Kutta methods: one explicit (ERK) and the other
diagonally implicit (DIRK).  Through appropriately partitioning the
ODE right-hand side into explicit and implicit components
:eq:`ARKODE_ODE_split_linearly_implicit`, such methods have the potential to
enable accurate and efficient time integration of stiff, nonstiff, and
mixed stiff/nonstiff systems of ordinary differential equations.  A
key feature allowing for high efficiency of these methods is that only
the components in :math:`f^I(t,y)` must be solved implicitly, allowing
for splittings tuned for use with optimal implicit solver algorithms.

This framework allows for significant freedom over the constitutive
methods used for each component, and ARKODE is packaged with a wide
array of built-in methods for use.  These built-in Butcher tables
include adaptive explicit methods of orders 2-9, adaptive implicit
methods of orders 2-5, and adaptive ImEx methods of orders 2-5.


*ERKStep* focuses specifically on problems posed in explicit form,

.. math::
   \dot{y} = f(t,y),  \qquad y(t_0) = y_0.
   :label: ARKODE_ODE_explicit

allowing for increased computational efficiency and memory savings.
The algorithms used in ERKStep are adaptive- and fixed-step explicit
Runge--Kutta methods.   As with ARKStep, the ERKStep module is packaged
with adaptive explicit methods of orders 2-9.

*SPRKStep* focuses on Hamiltonian systems posed in the form,

.. math::
   H(t, p, q) = T(t, p) + V(t, q)

.. math::
   \dot{p} = f_1(t,q) = \frac{\partial V(t,q)}{\partial q}, \quad
   \dot{q} = f_2(t,p) = \frac{\partial T(t,p)}{\partial p},
   :label: ARKODE_ODE_hamiltonian

allowing for conservation of quadratic invariants.

*MRIStep* focuses specifically on problems posed in additive form,

.. math::
   \dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y), \qquad y(t_0) = y_0.
   :label: ARKODE_ODE_two_rate

where here the right-hand side function is additively split into three
components:

* :math:`f^E(t,y)` contains the "slow-nonstiff" components of the system
  (this will be integrated using an explicit method and a large time step
  :math:`h^S`),

* :math:`f^I(t,y)` contains the "slow-stiff" components of the system
  (this will be integrated using an implicit method and a large time step
  :math:`h^S`), and

* :math:`f^F(t,y)` contains the "fast" components of the system (this will be
  integrated using a possibly different method than the slow time scale and a
  small time step :math:`h^F \ll h^S`).

For such problems, MRIStep provides fixed-step slow step multirate infinitesimal
step (MIS), multirate infinitesimal GARK (MRI-GARK), and implicit-explicit
MRI-GARK (IMEX-MRI-GARK) methods, allowing for evolution of the problem
:eq:`ARKODE_ODE_two_rate` using multirate methods having orders of accuracy 2-4.

For ARKStep or MRIStep problems that include nonzero implicit term
:math:`f^I(t,y)`, the resulting implicit system (assumed nonlinear, unless
specified otherwise) is solved approximately at each integration step, using a
SUNNonlinearSolver module, supplied either by the user or from the underlying
SUNDIALS infrastructure.  For nonlinear solver algorithms that internally
require a linear solver, ARKODE may use a variety of SUNLinearSolver modules
provided with SUNDIALS, or again may utilize a user-supplied module.


Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: ../../../shared/RecentChanges.rst

For changes in prior versions of SUNDIALS see :numref:`Changelog`.


Reading this User Guide
=======================

This user guide is a combination of general usage instructions and
specific example programs.  We expect that some readers will want to
concentrate on the general instructions, while others will refer
mostly to the examples, and the organization is intended to
accommodate both styles.

The structure of this document is as follows:

* In the next section we provide a thorough presentation of the
  underlying :ref:`mathematical algorithms <ARKODE.Mathematics>` used within
  the ARKODE family of solvers.

* We follow this with an overview of how the source code for
  ARKODE is :ref:`organized <ARKODE.Organization>`.

* The largest section follows, providing a full account of how to use
  ARKODE within C and C++ applications, including any instructions that are
  specific to a given time-stepping modules, :ref:`ARKStep <ARKODE.Usage.ARKStep>`,
  :ref:`ERKStep <ARKODE.Usage.ERKStep>`, or :ref:`MRIStep <ARKODE.Usage.MRIStep>`.
  This section then includes additional information on how to use ARKODE from
  applications written in :ref:`Fortran <Fortran>`, as well as information
  on how to leverage :ref:`GPU accelerators within ARKODE <SUNDIALS.GPU>`.

* A much smaller section follows, describing ARKODE's
  :ref:`Butcher table structure <ARKodeButcherTable>`, that is used by
  both ARKStep and ERKStep.

* Subsequent sections discuss shared SUNDIALS features that are used
  by ARKODE:
  :ref:`vector data structures <NVectors>`,
  :ref:`matrix data structures <SUNMatrix>`,
  :ref:`linear solver data structures <SUNLinSol>`,
  :ref:`nonlinear solver data structures <SUNNonlinSol>`,
  :ref:`memory management utilities <SUNMemory>`,
  and the
  :ref:`installation procedure <Installation>`.

* The final sections catalog the full set of :ref:`ARKODE constants
  <ARKODE.Constants>`, that are used for both input specifications and return
  codes, and the full set of :ref:`Butcher tables <Butcher>` that are
  packaged with ARKODE.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.
