.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
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

.. _SUNStepper:

#####################################
Stepper Data Structure
#####################################

This section presents the :c:type:`SUNStepper` base class which represents a
generic solution procedure for IVPs of the form

.. math::
   \dot{v}(t) = f(t, v) + r(t), \qquad v(t_0) = v_0,
   :label: SUNStepper_IVP

on an interval :math:`t \in [t_0, t_f]`. The time dependent forcing term,
:math:`r_i(t)`, is given by

.. math::
   r(t) = \sum_{k = 0}^{n_{\text{forcing}}-1}
   \left( \frac{t - t_{\text{shift}}}{t_{\text{scale}}} \right)^{k} \widehat{f}_k.
   :label: SUNStepper_forcing

:c:type:`SUNStepper` provides an abstraction over SUNDIALS integrators, custom
integrators, exact solution procedures, or other approaches for solving
:eq:`SUNStepper_IVP`. These are used, for example, in operator splitting and
forcing methods to solve inner IVPs in a flexible way.
