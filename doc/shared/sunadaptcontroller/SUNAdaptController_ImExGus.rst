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

.. _SUNAdaptController.ImExGus:

The SUNAdaptController_ImExGus Module
======================================

The ImEx Gustafsson implementation of the SUNAdaptController class, SUNAdaptController_ImExGus,
implements a combination of two adaptivity controllers proposed
by K. Gustafsson.  His "explicit" controller was proposed in :cite:p:`Gust:91`,
is primarily useful with explicit Runge--Kutta methods, and has the form

.. math::
   h' \;=\; \begin{cases}
      h_1\; \varepsilon_1^{-1/(p+1)}, &\quad\text{on the first step}, \\
      h_n\; \varepsilon_n^{-k_1^E/(p+1)}\;
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{k_2^E/(p+1)}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: expGusController

Similarly, Gustafsson's "implicit" controller was proposed in :cite:p:`Gust:94`
with the form

.. math::
   h' = \begin{cases}
      h_1 \varepsilon_1^{-1/(p+1)}, &\quad\text{on the first step}, \\
      h_n \left(\dfrac{h_n}{h_{n-1}}\right) \varepsilon_n^{-k_1^I/(p+1)}
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{-k_2^I/(p+1)}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: impGusController

In the above formulas, the default values of :math:`k_1^E`, :math:`k_2^E`,
:math:`k_1^I`, and :math:`k_2^I` are 0.367, 0.268, 0.98, and 0.95, respectively,
and :math:`p` is the global order of the time integration method.

The SUNAdaptController_ImExGus controller implements both formulas
:eq:`expGusController` and :eq:`impGusController`, and sets its recommended step
size as the minimum of these two.  It is implemented as a derived SUNAdaptController
class, and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_ImExGus {
     sunrealtype k1e;
     sunrealtype k2e;
     sunrealtype k1i;
     sunrealtype k2i;
     sunrealtype bias;
     sunrealtype ep;
     sunrealtype hp;
     sunbooleantype firststep;
   };

These entries of the *content* field contain the following information:

* ``k1e, k2e`` - explicit controller parameters used in :eq:`expGusController`.

* ``k1i, k2i`` - implicit controller parameters used in :eq:`impGusController`.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``hp`` - storage for the previous step size, :math:`h_{n-1}`.

* ``firststep`` - flag indicating whether a step has completed successfully, allowing
  the formulas above to transition between :math:`h_1` and :math:`h_n`.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_imexgus.h``.


The SUNAdaptController_ImExGus class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_H`` controller listed in
:numref:`SUNAdaptController.Description.operations`. The
SUNAdaptController_ImExGus class also provides the following additional user-callable
routines:


.. c:function:: SUNAdaptController SUNAdaptController_ImExGus(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_ImExGus
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object; otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_ImExGus(sunctx);

.. c:function:: SUNErrCode SUNAdaptController_SetParams_ImExGus(SUNAdaptController C, sunrealtype k1e, sunrealtype k2e, sunrealtype k1i, sunrealtype k2i)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_ImExGus object.
   :param k1e: parameter used within the controller time step estimate.
   :param k2e: parameter used within the controller time step estimate.
   :param k1i: parameter used within the controller time step estimate.
   :param k2i: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_ImExGus(C, 0.4, 0.3, -1.0, 1.0);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_imexgus".
