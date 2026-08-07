.. ----------------------------------------------------------------
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

.. _SUNAdjoint.Introduction:

Introduction to Adjoint Sensitivity Analysis
============================================

This section presents the :c:type:`SUNAdjointStepper` and
:c:type:`SUNAdjointCheckpointScheme` classes. The :c:type:`SUNAdjointStepper`
represents a generic adjoint sensitivity analysis (ASA) procedure to obtain the adjoint
sensitivities of an IVP of the form

.. math::
   \dot{y}(t) = f(t, y, p), \qquad y(t_0) = y_0(p), \qquad y \in \mathbb{R}^N,
   :label: SUNADJOINT_IVP

where :math:`p` is some set of :math:`N_s` problem parameters.

.. note::
   The API itself does not implement ASA, but it provides a common
   interface for ASA capabilities implemented in the SUNDIALS packages. Right now it supports :ref:`the
   ASA capabilities in ARKODE <ARKODE.Mathematics.ASA>`, while the ASA capabilities in :ref:`CVODES
   <CVODES.Mathematics.ASA>` and :ref:`IDAS <IDAS.Mathematics.ASA>` must be used directly.

Suppose we have a functional :math:`g(t_f, y(t_f), p)` for which we would like to compute the gradients
:math:`dg(t_f, y(t_f), p)/dy(t_0)` and/or :math:`dg(t_f, y(t_f), p)/dp`.
This most often arises in the form of an optimization problem such as

.. math::
   \min_{y(t_0), p} g(t_f, y(t_f), p)
   :label: SUNADJOINT_OPTIMIZATION_PROBLEM


.. warning::
   The CVODES documentation uses :math:`\lambda` to represent the adjoint variables needed
   to obtain the gradient :math:`dG/dp` where :math:`G` is an integral of :math:`g`.
   Our use of :math:`\lambda` in the following is akin to the use of :math:`\mu` in the CVODES docs.

The adjoint method is one approach to obtaining the gradients that is particularly efficient when
there are relatively few functionals and a large number of parameters. While :ref:`CVODES
<CVODES.Mathematics.ASA>` and :ref:`IDAS <IDAS.Mathematics.ASA>` *continuous* adjoint methods
(differentiate-then-discretize), ARKODE provides *discrete* adjoint methods
(discretize-then-differentiate). For the continuous approach, we derive and solve the adjoint IVP
backwards in time

.. math::
   \dot{\lambda}(t) = -f_y^*(t, y, p) \lambda,\quad \lambda(t_F) = g_y^*(t_f, y(t_f), p)
   :label: SUNADJOINT_CONTINUOUS_ADJOINT_IVP

where :math:`\lambda(t) \in \mathbb{R}^{N_s}`,
:math:`f_y \equiv \partial f/\partial y \in \mathbb{R}^{N \times N}` and
:math:`g_y \equiv \partial g/\partial y \in \mathbb{R}^{N \times N}`,
are the Jacobians with respect to the dependent variable, :math:`*` denotes the
Hermitian (conjugate) transpose, :math:`N` is the size of the original IVP, and
:math:`N_s` is the number of parameters. When solved with a numerical time
integration scheme, the solution to the continuous adjoint IVP is a numerical
approximation of the continuous adjoint sensitivities,

.. math::
   \lambda(t_n) \approx g_y(t_f, y(t_n), p), \quad \lambda(t_0) \approx g_y(t_f, y(t_0), p).
   :label: SUNADJOINT_CONTINUOUS_ADJOINT_SOLUTION

The gradients with respect to the parameters can then be obtained as

.. math::
   \frac{d g(t_f, y(t_n), p)}{dp} = \lambda^*(t_n) y_p(t_n) + g_p(t_f, y(t_n), p) + \int_{t_n}^{t_f} \lambda^*(t) f_p(t, y(t_n), p)~ dt,
   :label: SUNADJOINT_CONTINUOUS_PARAMETER_GRADIENT

where :math:`y_p(t) \equiv \partial y(t)/\partial p \in \mathbb{R}^{N \times N_s}`, and
:math:`g_p \equiv \partial g/\partial p \in \mathbb{R}^{N \times N_s}` and
:math:`f_p \equiv \partial f/\partial p \in \mathbb{R}^{N \times N_s}` are the
Jacobians with respect to the parameters.

For the discrete adjoint approach, we first numerically discretize the original IVP :eq:`SUNADJOINT_IVP`
using a time integration scheme, :math:`\varphi`, so that

.. math::
   y_0 = y(t_0),\quad y_n = \varphi(y_{n-k}, \cdots, y_{n-1}, p), \quad k = n, \cdots, 1.
   :label: SUNADJOINT_DISCRETE_IVP

For linear multistep methods :math:`k \geq 1` and for one step methods :math:`k = 1`.
Reformulating the optimization problem for the discrete case, we have

.. math::
   \min_{y_0, p} g(t_f, y_n, p)
   :label: SUNADJOINT_DISCRETE_OPTIMIZATION_PROBLEM

The gradients of :eq:`SUNADJOINT_DISCRETE_OPTIMIZATION_PROBLEM` can be computed using the transposed chain
rule backwards in time to obtain the discrete adjoint variables :math:`\lambda_n, \lambda_{n-1}, \cdots, \lambda_0`
and :math:`\mu_n, \mu_{n-1}, \cdots, \mu_0`.
The discrete adjoint variables represent the gradients of the discrete cost function
:eq:`SUNADJOINT_DISCRETE_OPTIMIZATION_PROBLEM` with respect to
changes in the discretized IVP :eq:`SUNADJOINT_DISCRETE_IVP`,

.. math::
   \frac{dg}{dy_n} = \lambda_n , \quad \frac{dg}{dp} = \mu_n + \lambda_n^* \left(\frac{\partial y_0}{\partial p} \right).
   :label: SUNADJOINT_DISCRETE_ADJOINT_GRADIENTS


.. _SUNAdjoint.DiscreteContinuous:

Discrete vs. Continuous Adjoint Method
--------------------------------------

It is understood that the continuous adjoint method can be problematic in the context of
optimization problems because the continuous adjoint method provides an approximation to the
gradient of a continuous cost function while the optimizer is expecting the gradient of the discrete
cost function. The discrepancy means that the optimizer can fail to due to inconsistent gradients
:cite:p:`giles2000introduction,gholami2019anode`. On the other hand, the discrete adjoint method
provides the exact gradient of the discrete cost function allowing the optimizer to fully converge.
Consequently, the discrete adjoint method is often preferable in optimization despite its own
drawbacks -- such as its (relatively) increased memory usage and the possible introduction of
unphysical computational modes :cite:p:`sirkes1997finite`. This is not to say that the discrete
adjoint approach is always the better choice over the continuous adjoint approach in optimization.
Computational efficiency and stability of one approach over the other can be both problem and method
dependent. Section 8 in the paper :cite:p:`rackauckas2020universal` discusses the tradeoffs further
and provides numerous references that may help inform users in choosing between the discrete and
continuous adjoint approaches.
