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

.. _SUNAdaptController.MRIHTol:

The SUNAdaptController_MRIHTol Module
======================================

.. versionadded:: 7.2.0

Mathematical motivation
-----------------------

The MRIHTol implementation of the SUNAdaptController class,
SUNAdaptController_MRIHTol, implements a general structure for telescopic
multirate temporal control.  A SUNAdaptController_MRIHTol object is constructed
using two single-rate controller objects, *HControl* and *TolControl*.  The
MRIHTol controller assumes that overall solution error at a given time scale
results from two types of error:

#. "Slow" temporal errors introduced at the current time scale,

   .. math::
      \varepsilon^S_{n} = C(t_n) \left(h_n^S\right)^{P+1},
      :label: slow_error_assumption

   where :math:`C(t)` is independent of the current time scale step size :math:`h_n^S`
   but may vary in time.

#. "Fast" errors introduced through calls to the next-fastest ("inner") solver,
   :math:`\varepsilon^F_{n}`.  If this inner solver is called to evolve IVPs over
   time intervals :math:`[t_{0,i}, t_{F,i}]` with a relative tolerance
   :math:`\text{RTOL}_n^F`, then it will result in accumulated errors over these
   intervals of the form

   .. math::
      \varepsilon^F_{n} = c(t_n) h_n^S \left(\text{RTOL}_n^F\right),

   where :math:`c(t)` is independent of the tolerance or subinterval width but may vary in
   time, or equivalently,

   .. math::
      \varepsilon^F_{n} = \kappa(t_n) \left(\text{tolfac}_n^F\right),
      :label: inner_solver_assumption

   where :math:`\text{RTOL}_n^F = \text{RTOL}^S \text{tolfac}_n^F`,
   the relative tolerance that was supplied to the current time scale
   solver is :math:`\text{RTOL}^S`, and
   :math:`\kappa(t_n) = c(t_n) h_n^S \text{RTOL}^S` is
   independent of the relative tolerance factor, :math:`\text{tolfac}_n^F`.

Single-rate controllers are constructed to adapt a single parameter, e.g.,
:math:`\delta`, under an assumption that solution error :math:`\varepsilon` depends
asymptotically on this parameter via the form

.. math::
   \varepsilon = \mathcal{O}(\delta^{q+1}).

Both :eq:`slow_error_assumption` and :eq:`inner_solver_assumption` fit this form,
with control parameters :math:`h_n^S` and :math:`\text{tolfac}^F_n`, and "orders"
:math:`P` and :math:`0`, respectively.  Thus an MRIHTol controller employs
*HControl* to adapt :math:`h_n^S` to control the current time scale error
:math:`\varepsilon^S_n`, and it employs *TolControl* to adapt
:math:`\text{tolfac}_n^F` to control the accumulated inner solver error
:math:`\varepsilon^F_n`.

To avoid overly large changes in calls to the inner solver, we apply bounds on the
results from *TolControl*.  If *TolControl* predicts a control parameter
:math:`\text{tolfac}'`, we obtain the eventual tolerance factor via
enforcing the following bounds:

.. math::
   \frac{\text{tolfac}_{n}^F}{\text{tolfac}'} &\le relch_{\text{max}},\\
   \frac{\text{tolfac}'}{\text{tolfac}_{n}^F} &\le relch_{\text{max}},\\
   \text{tolfac}_{min} &\le \text{tolfac}' \le \text{tolfac}_{max}.

The default values for these bounds are :math:`relch_{\text{max}} = 20`,
:math:`\text{tolfac}_{min} = 10^{-5}`, and :math:`\text{tolfac}_{max} = 1`.


Implementation
--------------

The SUNAdaptController_MRIHTol controller is implemented as a derived
:c:type:`SUNAdaptController` class, and its *content* field is defined by
the :c:struct:`SUNAdaptControllerContent_MRIHTol_` structure:

.. c:struct:: SUNAdaptControllerContent_MRIHTol_

   The member data structure for an MRIHTol controller

   .. c:member:: SUNAdaptController HControl

      A single time-scale controller to adapt the current step size, :math:`h^S_n`.

   .. c:member:: SUNAdaptController TolControl

      A single time-scale controller to adapt the inner solver relative tolerance
      factor, :math:`\text{reltol}^F_n`.

   .. c:member:: sunrealtype inner_max_relch

      The parameter :math:`relch_{\text{max}}` above.

   .. c:member:: sunrealtype inner_min_tolfac

      The parameter :math:`\text{tolfac}_{min}` above.

   .. c:member:: sunrealtype inner_max_tolfac

      The parameter :math:`\text{tolfac}_{max}` above.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_mrihtol.h``.

The SUNAdaptController_MRIHTol class provides implementations of all operations
relevant to a :c:enumerator:`SUN_ADAPTCONTROLLER_MRI_H_TOL` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_MRIHTol(SUNAdaptController HControl, SUNAdaptController TolControl, SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_MRIHTol
   object, and inserts its default parameters.

   :param HControl: the slow time step adaptivity controller object.
   :param TolControl: the inner solver tolerance factor adaptivity controller object.
   :param sunctx: the current :c:type:`SUNContext` object.

   :returns: if successful, a usable :c:type:`SUNAdaptController` object;
             otherwise it will return ``NULL``.


.. c:function:: SUNErrCode SUNAdaptController_SetParams_MRIHTol(SUNAdaptController C, sunrealtype inner_max_relch, sunrealtype inner_min_tolfac, sunrealtype inner_max_tolfac)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.  If any argument is outside the allowable range, that parameter
   will be reset to its default value.

   :param C: the SUNAdaptController_MRIHTol object.
   :param inner_max_relch: the parameter :math:`relch_{\text{max}}` (must be :math:`\ge 1`).
   :param inner_min_tolfac: the parameter :math:`\text{tolfac}_{min}` (must be :math:`> 0`).
   :param inner_max_tolfac: the parameter :math:`\text{tolfac}_{max}` (must be :math:`> 0`).

   :returns: :c:type:`SUNErrCode` indicating success or failure.
   
   .. versionchanged:: 7.5.0
   
      Removed the requirement that ``inner_max_tolfac`` must be :math:`\le 1`
   

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_mrihtol".


Usage
-----

Since this adaptivity controller is constructed using multiple single-rate adaptivity
controllers, there are a few steps required when setting this up in an application
(the steps below in *italics* correspond to the surrounding steps described in the
:ref:`MRIStep usage skeleton <ARKODE.Usage.MRIStep.Skeleton>`.

#. *Create an inner stepper object to solve the fast (inner) IVP*

#. Configure the inner stepper to use temporal adaptivity.  For example, when using
   an ARKODE inner stepper and the :c:func:`ARKodeCreateMRIStepInnerStepper`
   function, then either use its default adaptivity approach or supply a
   single-rate SUNAdaptController object, e.g.

   .. code:: C

      void* inner_arkode_mem = ERKStepCreate(f_f, T0, y, sunctx);
      MRIStepInnerStepper inner_stepper = nullptr;
      retval = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
      SUNAdaptController fcontrol = SUNAdaptController_PID(sunctx);
      retval = ARKodeSetAdaptController(inner_arkode_mem, fcontrol);

#. If using an ARKODE inner stepper, then set the desired temporal error accumulation
   estimation strategy via a call to :c:func:`ARKodeSetAccumulatedErrorType`, e.g.,

   .. code:: C

      retval = ARKodeSetAccumulatedErrorType(inner_arkode_mem, ARK_ACCUMERROR_MAX);

#. *Create an MRIStep object for the slow (outer) integration*

#. Create single-rate controllers for both the slow step size and inner solver
   tolerance, e.g.,

   .. code:: C

      SUNAdaptController scontrol_H   = SUNAdaptController_PI(sunctx);
      SUNAdaptController scontrol_Tol = SUNAdaptController_I(sunctx);

#. Create the multirate controller object, e.g.,

   .. code:: C

      SUNAdaptController scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);

#. Attach the multirate controller object to MRIStep, e.g.,

   .. code:: C

      retval = ARKodeSetAdaptController(arkode_mem, scontrol);

An example showing the above steps is provided in
``examples/arkode/CXX_serial/ark_kpr_nestedmri.cpp``, where multirate controller objects
are used for both the slow and intermediate time scales in a 3-time-scale simulation.
