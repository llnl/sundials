.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
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

.. _ARKODE.Usage.MRIStep.CustomInnerStepper:

MRIStep Custom Inner Steppers
=============================

Recall that infinitesimal multirate methods require solving a set of auxiliary IVPs

.. math::
   \dot{v}(t) = f^F(t, v) + r_i(t), \qquad v(t_{i,0}) = v_{i,0},
   :label: ARKODE_MRI_IVP

on intervals :math:`t \in [t_{i,0}, t_{i,f}]`.  For the MIS, MRI-GARK and IMEX-MRI-GARK
methods implemented in MRIStep, the forcing term :math:`r_i(t)`
presented in :numref:`ARKODE.Mathematics.MRIStep` can be equivalently written as

.. math::
   r_i(t) =
   \sum\limits_{k \geq 1} \hat{\omega}_{i,k} \tau^{k-1}
   +
   \sum\limits_{k \geq 1} \hat{\gamma}_{i,k} \tau^{k-1}
   :label: ARKODE_MRI_forcing_poly

where :math:`\tau = (t - t_{n,i-1}^S)/(h^S \Delta c_i^S)` is the normalized time
with :math:`\Delta c_i^S=\left(c^S_i - c^S_{i-1}\right)`, the slow stage times are
:math:`t_{n,i-1}^S = t_{n-1} + c_{i-1}^S h^S`, and the polynomial coefficient
vectors are

.. math::
   \hat{\omega}_{i,k} = \frac{1}{\Delta c_i^S} \sum\limits_{j=1}^{i-1}
   \Omega_{i,j,k} f^E(t_{n,j}^S, z_j)
   \quad\text{and}\quad
   \hat{\gamma}_{i,k} = \frac{1}{\Delta c_i^S} \sum\limits_{j=1}^i
   \Gamma_{i,j,k} f^I(t_{n,j}^S, z_j).
   :label: ARKODE_MRI_forcing_coefficients

The MERK and IMEX-MRI-SR methods included in MRIStep compute the forcing polynomial
:eq:`ARKODE_MRI_forcing_poly` similarly, with appropriate modifications to
:math:`\Delta c_i^S`, :math:`t_{n,i-1}^S`, and the coefficients
:eq:`ARKODE_MRI_forcing_coefficients`.

To evolve the IVP :eq:`ARKODE_MRI_IVP` MRIStep utilizes a generic time integrator
interface defined by the :c:type:`MRIStepInnerStepper` base class. This section
presents the :c:type:`MRIStepInnerStepper` base class and methods that define
the integrator interface as well as detailing the steps for creating an
:c:type:`MRIStepInnerStepper`.

.. toctree::
   :maxdepth: 1

   Description
   Implementing
